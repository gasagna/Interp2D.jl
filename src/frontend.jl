# Linear barycentric interpolator
type Linear2DInterpolator
	x::Vector{Float64}
	y::Vector{Float64}
	z::Vector{Float64}
	p::Vector{Int64}
	iadj::Vector{Int64}
	iend::Vector{Int64}
	ier::Vector{Int64}      # Preallocate these three guys 
	ist::Vector{Int64}      # to avoid repeated allocation
	ztmp::Vector{Float64}  # 
	function Linear2DInterpolator(x::Vector{Float64}, 
								  y::Vector{Float64}, 
								  z::Vector{Float64}; 
								  copy::Bool=true)
		if copy == true
			x, y, z = x[:], y[:], z[:]
		end
		x, y, z, p = reordr!(x, y, z)
		iadj, iend = trmesh(x, y)
		new(x, y, z, p, iadj, iend, [3], [1], [0.0])
	end
end

# pointwise evaluation. This does not allocate memory and should be fast.
evaluate(linint::Linear2DInterpolator, px::Float64, py::Float64	) = 
		intrc0!(px, py, linint.ztmp, linint.x, linint.y, linint.z, linint.iadj, 
				linint.iend, linint.ier, linint.ist)

# If znew is provided, then it is used to interpolate data, instead of using the original z data
evaluate!(linint::Linear2DInterpolator, px::Float64, py::Float64, znew::Vector; permute::Bool=false) =
	      intrc0!(px, py, linint.ztmp, linint.x, linint.y, 
	      	      permute == true ? permute!(znew, linint.p) : znew, linint.iadj, 
	      	      linint.iend, linint.ier, linint.ist)

# evaluation for a vector of points. In place version.
evaluate!(linint::Linear2DInterpolator, xi::Vector, yi::Vector, zi::Vector) =
		  intrc0!(xi, yi, zi, linint.ztmp, linint.x, linint.y, linint.z, linint.iadj,
		          linint.iend, linint.ier, linint.ist)

# evaluation for a vector of points. Allocates memory.
evaluate(linint::Linear2DInterpolator, xi::Vector, yi::Vector) =
		 evaluate!(linint, xi, yi, similar(xi))

# If znew is provided, then it is used to interpolate data, instead of using the original z data
evaluate!(linint::Linear2DInterpolator, xi::Vector, yi::Vector, znew::Vector,  zi::Vector; permute::Bool=false) =
	      intrc0!(xi, yi, zi, linint.ztmp, linint.x, linint.y, 
		          permute == true ? permute!(znew, linint.p) : znew, linint.iadj, 
		          linint.iend, linint.ier, linint.ist)



# Cubic Clough-Tocher interpolator
type Cubic2DInterpolator
	x::Vector{Float64}
	y::Vector{Float64}
	z::Vector{Float64}
	p::Vector{Int64}
	iadj::Vector{Int64}
	iend::Vector{Int64}
	zxzy::Matrix{Float64}
	ier::Vector{Int64}      # Preallocate these three guys 
	ist::Vector{Int64}      # to avoid repeated allocation
	ztmp::Vector{Float64}  # 
	function Cubic2DInterpolator(x::Vector{Float64}, 
								 y::Vector{Float64}, 
								 z::Vector{Float64}; 
								 copy::Bool=true,
								 nit::Int64=400, 
								 eps::Float64=1e-6)
		if copy == true
			x, y, z = x[:], y[:], z[:]
		end
		x, y, z, p = reordr!(x, y, z)
		iadj, iend = trmesh(x, y)
		zxzy = gradg(x, y, z, iadj, iend, nit, eps) 
		new(x, y, z, p, iadj, iend, zxzy, [3], [1], [0.0])
	end
end

# pointwise evaluation. This does not allocate memory and should be fast.
evaluate(cubint::Cubic2DInterpolator, px::Float64, py::Float64	) = 
		intrc1!(px, py, cubint.ztmp, cubint.x, cubint.y, cubint.z, cubint.iadj, 
				cubint.iend, cubint.ier, cubint.ist, cubint.zxzy)

# If znew is provided, then it is used to interpolate data, instead of using the original z data
evaluate!(cubint::Cubic2DInterpolator, px::Float64, py::Float64, znew::Vector; permute::Bool=false) =
		  intrc1!(px, py, cubint.ztmp, cubint.x, cubint.y, 
		          permute == true ? permute!(znew, cubint.p) : znew, cubint.iadj, 
		          cubint.iend, cubint.ier, cubint.ist, cubint.zxzy)

# evaluation for a vector of points. In place version.
evaluate!(cubint::Cubic2DInterpolator, xi::Vector, yi::Vector, zi::Vector) =
		  intrc1!(xi, yi, zi, cubint.ztmp, cubint.x, cubint.y, cubint.z, cubint.iadj,
		          cubint.iend, cubint.ier, cubint.ist, cubint.zxzy)

# evaluation for a vector of points. Allocates memory.
evaluate(cubint::Cubic2DInterpolator, xi::Vector, yi::Vector) =
		 evaluate!(cubint, xi, yi, similar(xi))

# there is no equivalent of passing a znew vector to the cubic because we 
# would have to recompute the gradients