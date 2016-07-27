using Base.Test
import Interp2D: reordr!, trmesh, gradg, intrc1!

# generate dummy data
srand(1)
N = 200
x = rand(N)
y = rand(N)

# simple function with known gradient
z = 2*x + 3*y

# ~~~ Test reordr ~~~
xold = x[:]
yold = y[:]
zold = z[:]
x, y, z, p = reordr!(x, y, z)
@test issorted(x)
@test xold[p] == x
@test yold[p] == y
@test zold[p] == z
@test_throws ErrorException reordr!(x[1:2], y, z)


# ~~~ Test trmesh ~~~
# Do not have tests on this function, but we call it to have a check
iadj, iend = trmesh(x, y)
@test_throws ErrorException trmesh(x[1:2], y)


# ~~~ Test gradg ~~~
zxzy = gradg(x, y, z, iadj, iend, 200, 1e-6)
zx = vec(zxzy[1, :])
zy = vec(zxzy[2, :])
@test norm(zx - 2.0)/norm(zx) < 1e-5
@test norm(zy - 3.0)/norm(zy) < 1e-5
@test_throws ErrorException gradg(x, y, z, iadj, iend, 1, -2.0)

