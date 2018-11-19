function fixedpointmap(f, x_0, tolerance, maxiter)
    D(f) = x -> ForwardDiff.derivative(f, x)
    f_prime = D(f)
    x_old = x_0
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - f(x_old)/f_prime(x_old)
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    return (value = x_old, normdiff = normdiff, iter = iter)
end

# define a map and parameters
f(x) = (x-1)^3

sol = fixedpointmap(f, 0, 1.0E-8, 1000) 
println("Fixed point = $(sol.value), and |x_t+1 - x_t| = $(sol.normdiff) in $(sol.iter) iterations")
