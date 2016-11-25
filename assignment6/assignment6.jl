
#= 
CPSC 302 Assignment 6
Nicholas Hu
=#

using Plots
pyplot()
using LaTeXStrings

# Tridiagonal system solver (Ax = b for A tridiagonal) without pivoting
#    
# Parameters: 4 vectors of size n or n-1 -- 
#             the subdiagonal (n-1), diagonal (n), superdiagonal (n-1), and 'b' vector (n)
#
# Returns:    1 vector of size n -- x = A^{-1} b

function trisolve{T<:Number}(dl::Vector{T}, d::Vector{T}, du::Vector{T}, b::Vector{T})::Vector{T}
    n = length(d)
    x = Vector{T}(n)
    
    # Row reduction
    for k = 1:n-1
        l = dl[k] / d[k]  # Multiplier
        d[k+1] = d[k+1] - l * du[k]
        b[k+1] = b[k+1] - l * b[k]
    end
    
    # Backward substitution
    x[n] = b[n] / d[n]
    for k = n-1:-1:1
        x[k] = (b[k] - du[k] * x[k+1]) / d[k]
    end
    
    return x
end;

n = 10
dl, d, du = Vector{Float64}(-1:-1:-(n-1)), Vector{Float64}(3:3:3n), Vector{Float64}(-2:-1:-n)

A, b = Tridiagonal(dl, d, du), vec(randn(n, 1))

x_julia = A \ b
x = trisolve(dl, d, du, b)

@printf("Norm of difference between solutions: %e", norm(x - x_julia))

g = t -> (pi / 2)^2 * sin((pi / 2) * t)

N = 100
h = (1 - 0) / N
t = linspace(0 + h, 1, N)
dl, d, du = [-vec(ones(N-2, 1)); -2], 2 * vec(ones(N)), -vec(ones(N-1, 1))

v = trisolve(dl / h^2, d / h^2, du / h^2, g(t))
u = sin((pi / 2) * t)

@printf("Infinity-norm of difference between solutions: %e", norm(v - u, Inf))

t = [0.0, 1.0, 2.0]
z = exp([0.1, 0.9, 2])

x_1, x_2 = linreg(t, log(z))
latexstring("u(t) = ", exp(x_1), "e^{", x_2, "t}")

m, n = 100, 10
A, b = randn(m, n), randn(m, 1)

lambda = linspace(0, 1, 50)
norm_res, norm_sol = Vector{Float64}(50), Vector{Float64}(50)

for i = 1:50
    x = [A; sqrt(lambda[i]) * eye(n)] \ [b; zeros(n, 1)]
    norm_res[i], norm_sol[i] = norm(b - A*x), norm(x)
end

p1 = plot(norm_sol, norm_res, label="", 
          xlabel=L"||\hat{x}||", 
          ylabel=L"||\vec{r}|| = ||\vec{b} - A\hat{x}||",
          title="Norm of residual vs. norm of solution")
p2 = plot(lambda, norm_sol, label="",
          xlabel=L"\lambda", ylabel=L"||\hat{x}||",
          title="Norm of solution vs. regularization parameter")
plot(p1, p2, layout=(2, 1), size=(800, 650))
