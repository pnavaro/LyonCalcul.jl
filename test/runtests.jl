using Test
using LyonCalcul

@testset "Test advection sinus function" begin

     xmin, xmax, nx = 0.0, 2π, 128
     dx = (xmax - xmin) / nx
     x = range(xmin, stop=xmax, length=nx+1)[1:end-1] |> collect
     f = sin.(x)
     advection! = Advection( nx, 5, dx )
     advection!( f, 0.5)
     @test maximum(abs.(f .- sin.(x .- 0.5))) ≈ 0.0 atol = 1e-12

end
