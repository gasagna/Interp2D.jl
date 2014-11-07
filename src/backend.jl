function reordr!(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
	length(x) == length(y) == length(z) || error("input vectors must be same size")
	N = int64(length(x))
	iflag = int64(3)    # reorder all vectors
	p = zeros(Int64, N) # permutation array
	ccall((:reordr_, lib624), 
	        Ptr{Void}, 
	        (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}), 
	        &N, &iflag, x, y, z, p)
	x, y, z, p
end

function trmesh(x::Vector{Float64}, y::Vector{Float64})
	length(x) == length(y) || error("input vectors must be same size")
	N = int64(length(x))
	iadj = zeros(Int64, 6*N-9)
	iend = zeros(Int64, N)
	ier = Int64[3]
	ccall((:trmesh_, lib624), 
	      Ptr{Void}, 
	      (Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}), 
	      &N, x, y, iadj, iend, ier)
	ier[1] == 1 && error("length of vectors less than 3")
	ier[1] == 2 && error("nodes are collinear")
	iadj, iend
end

gradg(x, y, z, iadj, iend, nit, eps) = 
	gradg!(x, y, z, iadj, iend, zeros(Float64, 2, length(x)), nit, eps)

function gradg!(x::Vector{Float64}, 
			    y::Vector{Float64}, 
			    z::Vector{Float64},
			    iadj::Vector{Int64},
			    iend::Vector{Int64},
			    zxzy::Matrix{Float64},
			    nit::Int64,
			    eps::Float64)
	length(x) == length(y) == length(z) || error("input vectors must be same size")
	size(zxzy) == (2, length(x)) || error("wrong dimension on input zxzy")
	N = int64(length(x))
	ier = Int64[3]  # init to unused code
	ccall((:gradg_, lib624), 
      	  Ptr{Void}, 
      	  (Ptr{Int64},   # N
      	   Ptr{Float64}, # x
      	   Ptr{Float64}, # y
      	   Ptr{Float64}, # z
      	   Ptr{Int64},   # iadj
      	   Ptr{Int64},   # iend
      	   Ptr{Float64}, # eps
      	   Ptr{Int64},   # nit
		   Ptr{Float64}, # zxzy
      	   Ptr{Int64}),  # ier
      	   &N, x, y, z, iadj, iend, &eps, &int64(nit), zxzy, ier)
	ier[1] == 1 && warn("desired gradient convergence not reached in $nit iterations")
	ier[1] == 2 && error("eps or nit lower than zero, or N lower than 3")
	return zxzy
end

gradl(x, y, z, iadj, iend) = gradl!(x, y, z, iadj, iend, zeros(Float64, 2, length(x)))

function gradl!(x::Vector{Float64}, 
			    y::Vector{Float64}, 
			    z::Vector{Float64},
			    iadj::Vector{Int64},
			    iend::Vector{Int64},
			    zxzy::Matrix{Float64})
	length(x) == length(y) == length(z) || error("input vectors must be same size")
	size(zxzy) == (2, length(x)) || error("wrong dimension on input zxzy")
	N = int64(length(x))
	ier = Int64[3]  # init to unused code
	dx = Float64[0.0]
	dy = Float64[0.0]
	for k = 1:N
		ccall((:gradl_, lib624), 
      	  	  Ptr{Void}, 
      	  	  (Ptr{Int64},   # N
      	  	   Ptr{Int64},   # k
      	   	   Ptr{Float64}, # x
      	   	   Ptr{Float64}, # y
      	   	   Ptr{Float64}, # z
      	   	   Ptr{Int64},   # iadj
      	   	   Ptr{Int64},   # iend
		   	   Ptr{Float64}, # dx
		   	   Ptr{Float64}, # dy
      	   	   Ptr{Int64}),  # ier
      	   &N, &k, x, y, z, iadj, iend, dx, dy, ier)
		zxzy[1, k] = dx[1]
		zxzy[2, k] = dy[1]
		ier[1] == -1 && warn("N or k out of range")
		ier[1] == -2 && error("nodes are collinear")
	end
	return zxzy
end

function intrc1!(px::Vector{Float64}, 
	             py::Vector{Float64}, 
	             pz::Vector{Float64}, 
			     x::Vector{Float64}, 
				 y::Vector{Float64}, 
				 z::Vector{Float64}, 
				 iadj::Vector{Int64}, 
				 iend::Vector{Int64}, 
				 zxzy::Matrix{Float64})
	length(x) == length(y) == length(z) || error("input vectors must be same size")
	size(zxzy) == (2, length(x)) || error("wrong dimension on input zxzy")
	N = int64(length(x))
	ier = Int64[3]       # initialise to unused code
	pztmp = Float64[0.0] # temporary array 
	iflag = int64(1)     # derivatives are provided
	ist = Int64[1]       # init to one, but it is updated at each call
	for i = 1:length(pz)
		ccall((:intrc1_, lib624), 
		      Ptr{Void}, 
	       	  (Ptr{Int64},  # N
	      	  Ptr{Float64}, # px
	      	  Ptr{Float64}, # py
	      	  Ptr{Float64}, # x
	      	  Ptr{Float64}, # y
	      	  Ptr{Float64}, # z
	          Ptr{Int64},   # iadj
	      	  Ptr{Int64},   # iend
	      	  Ptr{Int64},   # iflag = 1, derivatives are provided
	      	  Ptr{Float64}, # zxzy
	      	  Ptr{Int64},   # ist
	      	  Ptr{Float64}, # pz
	      	  Ptr{Int64}),  # ier
	          &N, &px[i], &py[i], x, y, z, iadj, iend, &iflag, zxzy, ist, pztmp, ier)
		#ier[1] ==  1 && warn("extrapolating out of convex hull of data")
		ier[1] == -1 && error("N, iflag orr ist out of range")
		ier[1] == -2 && error("nodes are collinear")
		pz[i] = pztmp[1]
	end
	return pz
end

function intrc0!(px::Vector{Float64}, 
	             py::Vector{Float64}, 
	             pz::Vector{Float64}, 
			     x::Vector{Float64}, 
				 y::Vector{Float64}, 
				 z::Vector{Float64}, 
				 iadj::Vector{Int64}, 
				 iend::Vector{Int64})
	length(x) == length(y) == length(z) || error("input vectors must be same size")
	N = int64(length(x)) 
	ier = Int64[3]       # initialise to unused code
	pztmp = Float64[0.0] # temporary array 
	ist = Int64[1]       # init to one, but it is updated at each call
	for i = 1:length(pz)
		ccall((:intrc0_, lib624), 
		      Ptr{Void}, 
	       	  (Ptr{Int64},  # N
	      	  Ptr{Float64}, # px
	      	  Ptr{Float64}, # py
	      	  Ptr{Float64}, # x
	      	  Ptr{Float64}, # y
	      	  Ptr{Float64}, # z
	          Ptr{Int64},   # iadj
	      	  Ptr{Int64},   # iend
	      	  Ptr{Int64},   # ist
	      	  Ptr{Float64}, # pz
	      	  Ptr{Int64}),  # ier
	          &N, &px[i], &py[i], x, y, z, iadj, iend, ist, pztmp, ier)
		#ier[1] ==  1 && warn("extrapolating out of convex hull of data")
		ier[1] == -1 && error("N, iflag orr ist out of range")
		ier[1] == -2 && error("nodes are collinear")
		pz[i] = pztmp[1]
	end
	return pz
end