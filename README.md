Interp2D
========

A package for interpolation of scattered data in two dimensions, using linear or cubic approximation. Tested on Julia 0.4.

Usage
-----
Suppose we have a scattered dataset in two-dimensions, as a set of three arrays `x`, `y` and `z`.
.  A linear interpolator is created with
```julia
lint = Linear2DInterpolator(x, y, z)
```
and interpolation at new points `xi`, `yi` is obtained as 
```julia
evaluate(lint, xi, yi)
```
or with 
```julia
evaluate!(lint, xi, yi, zi)
```
if `zi` has already been allocated. 

For cubic interpolation the interpolator object is created similarly, but two additional options are provided:
```julia
cint = Cubic2DInterpolator(x, y, z; nit=400, eps=1e-6)
```
The two options control the calculation of the gradient at the nodes of the triangulation, and are set to good default values.  Might want to try different combinations if needed. Interpolation on new data points is performed exactly in the same way as for the linear case.

Algorithm description
---------------------
This package wraps a subset of the routines from the Algorithm 624 package available in the toms directory on netlib. For both the linear and the cubic approaches, a Thiessen triangulation of the input data is first generated. For the linear case, on a new input point `(xi, yi)` where interpolation is required, the triangle containing the point is found and linear barycentric interpolation with the three nodes is performed. For the cubic interpolation, after the triangulation step, the gradient of the data at the nodes is calculated using a global method, (see references). For new input points, a cubic interpolation is performed using the values at the nodes of the triangle containing the point as well as the gradients.