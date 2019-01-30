var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#LyonCalcul.Advection",
    "page": "Documentation",
    "title": "LyonCalcul.Advection",
    "category": "type",
    "text": "Advection(n, p, delta)\n\nn :: Number of points.\np :: Spline degree.\ndelta :: space size step.\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.Advection-Tuple{Array{Float64,1},Float64}",
    "page": "Documentation",
    "title": "LyonCalcul.Advection",
    "category": "method",
    "text": "advection! = Advection( n, p, delta)\nadvection!( f, alpha )\n\nCreate a function to compute the interpolating spline \nof degree p of odd\ndegree of a 1D function f on a periodic uniform mesh, at\nall points x after a displacement alpha. \nInput f type is Vector{Float64} and is updated inplace.\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.bspline-Tuple{Int64,Int64,Float64}",
    "page": "Documentation",
    "title": "LyonCalcul.bspline",
    "category": "method",
    "text": "bspline(p, j, x)\n\nReturn the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j. Implemented recursively using the De Boor\'s Algorithm\n\nB_i0(x) = left\nbeginmatrix\n1  mathrmif  quad t_i  x  t_i+1 \n0  mathrmotherwise \nendmatrix\nright\n\nB_ip(x) = fracx - t_it_i+p - t_i B_ip-1(x) \n+ fract_i+p+1 - xt_i+p+1 - t_i+1 B_i+1p-1(x)\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.jl-1",
    "page": "Documentation",
    "title": "LyonCalcul.jl",
    "category": "section",
    "text": "Documentation for LyonCalcul.jlModules = [LyonCalcul]\nOrder   = [:type, :function]"
},

]}
