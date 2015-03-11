# Linear barycentric interpolator
type Linear2DInterpolator
	x::Vector{Float64}
	y::Vector{Float64}
	z::Vector{Float64}
	p::Vector{Int64}
	iadj::Vector{Int64}
	iend::Vector{Int64}
	function Linear2DInterpolator(x::Vector{Float64}, 
								  y::Vector{Float64}, 
								  z::Vector{Float64}; 
								  copy::Bool=true)
		if copy == true
			x, y, z = x[:], y[:], z[:]
		end
		x, y, z, p = reordr!(x, y, z)
		iadj, iend = trmesh(x, y)
		new(x, y, z, p, iadj, iend)
	end
end

evaluate(linint::Linear2DInterpolator, xi::Vector, yi::Vector) =
	evaluate!(linint, xi, yi, similar(xi))

evaluate!(linint::Linear2DInterpolator, xi::Vector, yi::Vector, zi::Vector) =
	intrc0!(xi, yi, zi, linint.x, linint.y, linint.z, linint.iadj, linint.iend)

# if z is provided, then it is used to interpolate data, instead of using the original z data
evaluate!(linint::Linear2DInterpolator, xi::Vector, yi::Vector, zi::Vector, znew::Vector) =
	intrc0!(xi, yi, zi, linint.x, linint.y, znew[linint.p], linint.iadj, linint.iend)

# Cubic Clough-Tocher interpolator
type Cubic2DInterpolator
	x::Vector{Float64}
	y::Vector{Float64}
	z::Vector{Float64}
	p::Vector{Int64}
	iadj::Vector{Int64}
	iend::Vector{Int64}
	zxzy::Matrix{Float64}
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
		new(x, y, z, iadj, iend, zxzy, p)
	end
end

evaluate(cubint::Cubic2DInterpolator, xi::Vector, yi::Vector) =
	evaluate!(cubint, xi, yi, similar(xi))

evaluate!(cubint::Cubic2DInterpolator, xi::Vector, yi::Vector, zi::Vector) =
	intrc1!(xi, yi, zi, cubint.x, cubint.y, cubint.z, cubint.iadj, cubint.iend, cubint.zxzy)

# if z is provided, then it is used to interpolate data, instead of using the original z data
evaluate!(cubint::Cubic2DInterpolator, xi::Vector, yi::Vector, zi::Vector, znew::Vector) =
	intrc1!(xi, yi, zi, cubint.x, cubint.y, cubint.z[cubint.p], cubint.iadj, cubint.iend, cubint.zxzy)