function newtonsmethod(p::Poly, x₀; tolerance = 1E-7,maxiter = 100)
    x_old = x₀
    normdiff = Inf
    iter = 1
    p_prime = polyder(p)
    while normdiff >= tolerance && iter <= maxiter
        x_new = x_old - p(x_old) / p_prime(x_old)
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter += 1
    end
    return x_old
end

p = Poly([2, -5, 2], :x)
x₀ = 0.0
@show newtonsmethod(p,x₀)
@show roots(p)
