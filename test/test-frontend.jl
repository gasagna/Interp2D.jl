using Base.Test
using Interp2D

# example data
N = 10000
x = rand(N)
y = rand(N)
z = x + y

# test linear interpolator
lint = Linear2DInterpolator(x, y, z)
cint = Cubic2DInterpolator(x, y, z)

for interpolator in [lint, cint]
	# test for single point
	@test_approx_eq_eps evaluate(interpolator, 0.5,  0.5)  1.0 1e-8

	# test for vector inputs
	out = evaluate(interpolator, [0.5, 0.2],  [0.5, 0.8])
	@test_approx_eq_eps out [1.0, 1.0] 1e-8

	# test for vector input. Inplace version
	zi = zeros(Float64, 2)
	zi = evaluate!(interpolator, [0.5, 0.2],  [0.5, 0.8], zi)
	@test_approx_eq_eps zi [1.0, 1.0] 1e-8

	# test for vector input. pass new data.
	# only for linear interpolator
	if interpolator == lint
		znew = 2x + 2y
		zi = zeros(Float64, 2)
		zi = evaluate!(interpolator, [0.5, 0.2],  [0.5, 0.8], znew, zi; permute=true) # this call is permuting zew
		@test_approx_eq_eps zi [2.0, 2.0] 1e-8
	end
end


