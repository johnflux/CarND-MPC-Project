# The Model

The model consists of:

**The state variables**:

* **`x`**, **`y`**:
  The position of the car, in car coordinates
* **`psi`**:
  The angle of the car, in car coordinates, in radians.  0 is facing in the x direction
* **`v`**:
  The speed of the car

Since `x`, `y` and `psi` are in car coordinates, they are 0,0,0  at any given moment in time.  But for future model planning they are taken from the current car coordinates.

**Error terms**:

* **`cte`**: The cross track error.  The distance from the desired position (from waypoints) and the position (planned or actual).
* **`epsi`**: The angle error.  The error in the desired angle and angle (planned or actual).

**Actuators**:

* **`delta`**: The turning angle, in radians, between -25 degrees to 25 degrees.
* **`a`**: The throttle (acceleration).

## Update Equations

The update model used was the bicycle model:

    x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * d
    y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
    psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
    v_[t] = v[t-1] + a[t-1] * dt
    cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
    epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt

# Timestep

My main restriction was that my laptop really struggled to run the Unity app.  The computational restriction meant that I had to keep the length and duration low, to minimize the time spent in the solver, and thus minimize the latency.

Particularly at high speeds, the latency of the solver was an absolute killer.

With that in mind, I did try to reduce the number of Timesteps (`N`) and/or increase the timestep duration (`dt`), but found that the gives values of:

* `N` = 25
* `dt` = 0.05 seconds
* (`Total time = N*dt = 1.25 seconds`)

Worked best.

# MPC Preprocessing

I record the amount of time spent in total in the function, including the solver time and the simulated latency, and use a running average (to provide a high pass filter) to give an estimate total latency.

I use that total latency as the '`dt`' value in the model update equations mentioned above to predict the position of the car by the time the actuators are actuated.

I transform the MPC waypoints into this latency-adjusted predicted car coordinate system.

# Polynomial Fitting

The MPC waypoints, in the predicted car coordinate system, are polyfitted with a polynomial of order 3.

The polynomial is evaluated at `x=0`  (i.e. `coeffs[1]`) and subtracted from the current `x` (which is `0` in car coordinates) to give the **`cte`** and its differential is evaluated at `x=0` and subtracted from `0`, likewise, to give the **`epsi`**.

# Model Predictive Control with Latency

I mentioned latency above, but just to recap:

I measure the actual total latency (with a high pass filter) and apply the model by that latency to predict the car state.  This is used to preprocess the waypoints.

