#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
const int N = 15 ;
const double dt = 0.1;

const double ref_v_turning = 90;
const double ref_v_straight = 200;

// State
const int x_start = 0*N;
const int y_start = 1*N;
const int psi_start = 2*N;
const int v_start = 3*N;
const int cte_start = 4*N;
const int epsi_start = 5*N;
// Actuators
const int delta_start = 6*N; /* Size N-1 */
const int a_start = 7*N - 1; /* Size N-1 */
const int n_vars = 8*N - 2;

// Store the last good solution, so that we can reuse it if we have a bad solution
CppAD::vector<double> last_good_solution;
int last_good_solution_age = -1;

AD<double> polyeval_cppad(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}

AD<double> polyeval_gradiant_cppad(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * CppAD::pow(x, i-1);
  }
  return result;
}

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    const double curving = coeffs[1]*coeffs[1] + coeffs[2]*coeffs[2] + coeffs[3]*coeffs[3];
    const bool is_turning = curving > 0.001;
    const double ref_v = is_turning?ref_v_turning:ref_v_straight;
    for (int t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += 0.1*CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += 100*CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += 100 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (int t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = polyeval_cppad(coeffs, x0);
      AD<double> psides0 = polyeval_gradiant_cppad(coeffs, x0);

      // The equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, vector<double> &predicted_trajectory_x, vector<double> &predicted_trajectory_y) {
  bool ok = true;

  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  const int state_size = 6;
  const int num_actuators = 2;
  assert(state.size() == state_size);

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  assert(n_vars == state_size * N + num_actuators * (N - 1));
  // TODO: Set the number of constraints
  const int n_constraints = N * state_size;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    // Reuse last solution, if good, shifting by 1 time step
    if (last_good_solution_age == 0 && i != n_vars-1)
      vars[i] = last_good_solution[i+1];
    else
      vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
    for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // TODO: Set lower and upper limits for variables.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332; /* -25 degrees, in radians */
    vars_upperbound[i] = 0.436332;  /*  25 degrees, in radians */
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }



  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  constraints_upperbound[x_start] = constraints_lowerbound[x_start] = x;
  constraints_upperbound[y_start] = constraints_lowerbound[y_start] = y;
  constraints_upperbound[psi_start] = constraints_lowerbound[psi_start] = psi;
  constraints_upperbound[v_start] = constraints_lowerbound[v_start] = v;
  constraints_upperbound[cte_start] = constraints_lowerbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = constraints_lowerbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  2\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  if(last_good_solution_age > 0) // If we didn't find a solution last time, let's spend a bit more time on it this time
    options += "Numeric max_cpu_time          0.5\n";
  else
    options += "Numeric max_cpu_time          0.3\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok = solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  int offset = 0;
  CppAD::vector<double> solution_x;
  if(!ok) {
    std::cout << "SOLUTION IS BAD!" << endl;
    assert(last_good_solution_age!= -1); // A bad solution on the very first frame
    solution_x = last_good_solution;
    last_good_solution_age++;
    offset = last_good_solution_age;
  } else {
    solution_x = solution.x;
    last_good_solution = solution_x;
    last_good_solution_age = 0;
  }
  assert(offset < N); //If we get too many bad solutions, there's not much we can do

  predicted_trajectory_x.clear();
  predicted_trajectory_y.clear();
  for(int i = offset; i < N; ++i) {
    predicted_trajectory_x.push_back(solution_x[x_start+i]);
    predicted_trajectory_y.push_back(solution_x[y_start+i]);
  }

  // Cost
  //auto cost = solution.obj_value;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {solution_x[delta_start+offset], solution_x[a_start+offset]};
}
