using Base.Test
using Interp2D

# example data
N = 100
x = rand(N)
y = rand(N)
z = x + y

# test linear interpolator
lint = Linear2DInterpolator(x, y, z)
@test_approx_eq_eps evaluate(lint, [0.5],  [0.5])  [1.0] 1e-8
zi = [0.0]
evaluate!(lint, [0.5],  [0.5], zi)
@test_approx_eq_eps zi [1.0] 1e-8


# test cubic interpolator
cint = Cubic2DInterpolator(x, y, z)
@test_approx_eq_eps evaluate(cint, [0.5],  [0.5])  [1.0] 1e-8
zi = [0.0]
evaluate!(cint, [0.5],  [0.5], zi)
@test_approx_eq_eps zi [1.0] 1e-8


