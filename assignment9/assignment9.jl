
#=
CPSC 302 Assignment 9
Nicholas Hu
=#

# Newton's method for nonlinear systems
#
# Parameters: a vector-valued vector function f,
#             the Jacobian of f as a matrix-valued vector function J
#             an initial solution guess x0
#             a tolerance tol
#             the maximum number of iterations to perform (default 100)

function newton{T<:Number}(f::Function, J::Function, x0::Vector{T};
                           tol::Float64=eps(), maxits::Int=100)
    x::Vector{T} = x0
    p = Vector{T}(length(x))
    
    for k = 1:maxits
        p = J(x) \ f(x)
        x -= p
        norm(p) < tol * (1 + norm(x)) && return x
    end
end;

function func(x::Vector{Float64})::Vector{Float64}
    return [2x[1] - x[2] - exp(-x[1]);
            -x[1] + 2x[2] - exp(-x[2])]
end
    
function jac(x::Vector{Float64})::Matrix{Float64}
    return [2 + exp(-x[1]) -1;
            -1 2 + exp(-x[2])]
end;

x = newton(func, jac, [-5., -5.])
@printf("x = [%.16f, %.16f]\n", x[1], x[2])
@printf("||f(x)|| = %e\n", norm(func(x)))

y = -5
for i = 1:10
    y = (1 + y) / (1 + exp(y))
end

println("y = $y")
