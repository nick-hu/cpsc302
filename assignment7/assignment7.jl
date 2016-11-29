
#= 
CPSC 302 Assignment 7
Nicholas Hu
=#

using Plots
pyplot()
using LaTeXStrings

# Linear stationary iterative (i.e., relaxation) method
#
# Parameters: an n-by-n system matrix A,
#             a right-hand side n-vector b,
#             an n-by-n partial system matrix M,
#             an n-vector x to iterate upon (in-place!), and
#             the number of iterations to perform (default 1)

function relax!{T<:Number}(A::Matrix{T}, b::Vector{T}, M::Matrix{T}, 
                           x::Vector{T}, n::Int=1)
    for i = 1:n
        r = b - A*x
        x[:] = x[:] + M \ r
    end
end

jacobi!{T<:Number}(A::Matrix{T}, b::Vector{T}, x::Vector{T}, 
                   n::Int=1) = relax!(A, b, diagm(diag(A)), x, n)
gauss_seidel!{T<:Number}(A::Matrix{T}, b::Vector{T}, x::Vector{T},
                         n::Int=1) = relax!(A, b, tril(A), x, n);

A = ones(3, 3) + 9 * eye(3)
b = 12 * ones(3)
sol = ones(3)
err = x -> norm(sol - x, 1)

x = zeros(3)

@printf("||x* - x_0||_1 = %e\n", err(x))
jacobi!(A, b, x)
@printf("||x* - x_1||_1 = %e\n", err(x))
jacobi!(A, b, x)
@printf("||x* - x_2||_1 = %e\n", err(x))
gauss_seidel!(A, b, x)
@printf("||x* - x_3||_1 = %e\n", err(x))

spradius = T -> maximum(abs(eigvals(T)))
T = eye(3) - inv(diagm(diag(A))) * A
@printf("Spectral radius of T: %f\n", spradius(T))

A = 5 * ones(3, 3) - 3 * eye(3)
x = zeros(3)

@printf("||b - Ax_0||_1 = %e\n", err(x))
jacobi!(A, b, x)
@printf("||b - Ax_1||_1 = %e\n", err(x))
jacobi!(A, b, x)
@printf("||b - Ax_2||_1 = %e\n", err(x))

T = eye(3) - inv(diagm(diag(A))) * A
@printf("Spectral radius of T: %f\n", spradius(T))

a = 0.55
A = a * ones(3, 3) + (1 - a) * eye(3)
b = randn(3)
err = x -> norm(b - A*x)

x_jac = randn(3)
x_GS = x_jac[:]
err_jac, err_GS = Vector{Float64}(51), Vector{Float64}(51)

for k in 1:51
    err_jac[k], err_GS[k] = err(x_jac), err(x_GS)
    jacobi!(A, b, x_jac)
    gauss_seidel!(A, b, x_GS)
end;

p1 = plot(0:50, err_jac, yscale=:log10, label="",
          xlabel=L"k",
          ylabel=L"||\vec{r}^{(k)}|| = ||\vec{b} - A\vec{x}^{(k)}||",
          title="Norm of error vs. iteration number (Jacobi)")
p2 = plot(0:50, err_GS, yscale=:log10, label="",
          xlabel=L"k",
          ylabel=L"||\vec{r}^{(k)}|| = ||\vec{b} - A\vec{x}^{(k)}||",
          title="Norm of error vs. iteration number (Gauss-Seidel)")
plot(p1, p2, layout=(2, 1), size=(800, 650))
