export Advection

"""
    Advection(n, p, delta)

    n :: Number of points.
    p :: Spline degree.
    delta :: space size step.
    
"""
mutable struct Advection
    n         :: Int64
    p         :: Int64
    delta     :: Float64
    modes     :: Vector{Float64}
    eig_bspl  :: Vector{ComplexF64}
    eig_alpha :: Vector{ComplexF64}
    
    function Advection( n, p, delta )
        
        eig_bspl  = zeros(ComplexF64,n)
        eig_alpha = similar(eig_bspl)
        modes = 2π * (0:n-1) / n
   
        eig_bspl .= bspline(p, -div(p+1,2), 0.0)
        for j in 1:div(p+1,2)-1
            eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)
             * 2 * cos.(j * modes))
        end
        
        new(n, p, delta, modes, eig_bspl, eig_alpha)
    end
end

using FFTW

"""
    advection! = Advection( n, p, delta)
    advection!( f, alpha )

    Create a function to compute the interpolating spline 
    of degree p of odd
    degree of a 1D function f on a periodic uniform mesh, at
    all points x after a displacement alpha. 
    Input f type is Vector{Float64} and is updated inplace.


"""
function (adv :: Advection)(f     :: Vector{Float64},
                            alpha :: Float64)

  
   ishift = floor(- alpha / adv.delta)
   beta = - ishift - alpha / adv.delta
   fill!(adv.eig_alpha, 0.0)
   for j in -div(adv.p-1,2):div(adv.p+1,2)
      adv.eig_alpha .+= (bspline(adv.p, j-div(adv.p+1,2), beta)
         .* exp.((ishift + j) * 1im .* adv.modes))
   end

   f .= real(ifft(fft(f) .* adv.eig_alpha ./ adv.eig_bspl))

end
