@unix_only begin
    cd(joinpath(dirname(@__FILE__), "src"))
    run(`gfortran -freal-4-real-8 -finteger-4-integer-8 -O3 --shared -fPIC -o lib624.so 624.f`)
end