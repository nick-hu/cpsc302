
#= 
CPSC 302 Assignment 5
Nicholas Hu
=#

using Plots
pyplot()
using LaTeXStrings

# Function-array 'hybrid' type
immutable LSFunction <: Function
    coeffs::Vector{Float64}
    f::Function
    
    # Performs a least-squares fit given x, y, and basis functions
    # By default, the basis functions are set such that the result is a linear regression
    function LSFunction{F<:Function}(x, y, basis::Array{F}=[x->1, x->x])
        m, n = length(x), length(basis)
        A = Matrix{Float64}(m, n)
        x = Vector{Float64}(x)

        for j = 1:n
            A[:, j] = basis[j](x)
        end

        c = A \ y
        new(c, x -> sum(c[j] * basis[j](x) for j = 1:length(c)))
    end
end

# Function-like behaviour
(obj::LSFunction)(x) = obj.f(x)
# Array-like behaviour
Base.getindex(obj::LSFunction, i) = obj.coeffs[i]
Base.endof(obj::LSFunction) = endof(obj.coeffs)

# Generates an array of monomial basis functions for use in LSFunction
monomial(n::Int)::Array{Function} = [x->x.^j for j = 0:n]

t = 0:0.1:1.3
b = [0.95, 1.01, 1.05, 0.97, 0.0, -0.1, 0.02, 
     -0.1, 0.01, -0.15, 0.72, 0.79, 0.91, 1.0]

scatter(t, b; label="", xlabel=L"t", ylabel=L"b", 
        title=L"Material property vs. $t$", legendfont=font(11))

# Breakpoints
piece1, piece2, piece3 = 1:4, 5:10, 11:length(t)

const1 = LSFunction(t[piece1], b[piece1], monomial(0))
const2 = LSFunction(t[piece2], b[piece2], monomial(0))
line = LSFunction(t[piece3], b[piece3])

format = s -> @sprintf("%.4f", s)
plot!(t[piece1], const1, label=latexstring("b = $(format(const1[1]))"))
plot!(t[piece2], const2, label=latexstring("b = $(format(const2[1]))"))
plot!(t[piece3], line, 
      label=latexstring("b = $(format(line[2]))t $(format(line[1]))"))

A, b = randn(20, 10), randn(20)

x_bs = A \ b
x_inv = inv(A'A)A'b
x_pinv = pinv(A)b

@printf("||x_bs - x_inv||_2   = %e\n", norm(x_bs - x_inv))
@printf("||x_bs - x_pinv||_2  = %e\n", norm(x_bs - x_pinv))
@printf("||x_inv - x_pinv||_2 = %e\n", norm(x_inv - x_pinv))

for test in 1:10
    delta = randn(10)
    delta /= norm(delta) ^ test
    if (norm(A*(x_bs + delta) - b)^2 <= norm(A*x_bs - b)^2)
        @printf("Test %d failed!\n", test)
    end
end

A, b = randn(100000, 100), randn(100000)
@time A \ b;
@time A \ b;

t = [1971, 1972, 1974, 1978, 1982, 1985, 1989, 
     1993, 1997, 1999, 2000, 2002, 2003]
N = [2250, 2500, 5000, 29000, 120000, 275000, 1180000, 3100000, 
     7500000, 24000000, 42000000, 220000000, 410000000]

p = scatter(t, N; label="", yscale=:log10,
            xlabel=latexstring(L"t", " (year)"), 
            ylabel=latexstring(L"N", " (number of transistors)"), 
            title="Number of transistors vs. year", 
            legendfont=font(11))

lfit = LSFunction(t-1970, log10(N))
N_pred = exp10(lfit(t-1970))

plot!(t, N_pred, label=latexstring(L"\log_{10} \, N = ", 
      "\$ $(format(lfit[1])) + $(format(lfit[2])) (t-1970) \$"))
display(p)
@printf("RMS Error in linear scale: %e\n", norm(log(N_pred) - log(N)) / sqrt(length(N)))
@printf("RMS Error in exponential scale: %e\n", norm(N_pred - N) / sqrt(length(N)))

@printf("Predicted number of transistors in 2015: %e\n", exp10(lfit(2015-1970)))

# Moore's Law with doubling every 1.5 years
plot!(t, N[1] * exp2((t - t[1]) / 1.5), label="Moore's Law (1.5 a. doubling)", line=:dash)
# Moore's Law with doubling every 2 years
plot!(t, N[1] * exp2((t - t[1]) / 2), label="Moore's Law (2 a. doubling)", line=:dash)
