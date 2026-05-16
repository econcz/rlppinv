# rlppinv 1.0.0

## Changes
* Synchronized the package with `rclsp` 1.0.0.
* Added support for the `cond_tolerance` argument, passed through to the
  underlying CLSP solver.

## Bug fixes
* Corrected upper- and lower-bound vector construction to match the canonical
  CLSP block order.
* Improved compatibility with NumPy-style bound ordering used by the Python
  implementation.

# rlppinv 0.3.0

## Bug fixes
* Updated CVXR integration for CVXR 1.8.x compatibility.

# rlppinv 0.2.0

## Bug fixes
* Corrected row-wise reconstruction of x from the solution vector z.
* Updated the minimum required version of **rclsp** to **>= 0.3.0**.
