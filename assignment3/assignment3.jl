
#= 
CPSC 302 Assignment 3
Nicholas Hu
=#

using Plots
pyplot()
using LaTeXStrings

    
function print_correctdigits(approx, exact)
    isnan(exact) || @printf("%30d", floor(-log10(abs(approx - exact))))
    @printf("\n")
end


function newton(f, df, x_0; atol=1e-10, ftol=NaN, max_iter = Inf, verbose=false, exact_val=NaN)
    x_n = x_0
    iteration = 0

    if verbose
        @printf("%5s%15s", "k", "x_k")
        isnan(ftol) || @printf("%15s", "f(x_k)")
        isnan(exact_val) || @printf("%30s", "Correct digits in frac. part")
        @printf("\n%5d%15f", iteration, x_n)
        isnan(ftol) || @printf("%15.2f", f(x_n))
        print_correctdigits(x_n, exact_val)
    end

    while iteration <= max_iter
        x_prev = x_n
        x_n = x_n - f(x_n) / df(x_n)
        iteration += 1

        if verbose
            @printf("%5d%15f", iteration, x_n)
            isnan(ftol) || @printf("%15.2e", f(x_n))
            print_correctdigits(x_n, exact_val)
        end

        ((abs(x_n - x_prev) < atol) || (abs(f(x_n)) < ftol)) && break
    end

    @printf("\n")
    return x_n
end
    
    
F = a -> (x -> x^3 - a)
D = x -> 3x^2

newton(F(0), D, 1, ftol=1e-8, verbose=true)
newton(F(2), D, 1, ftol=1e-8, verbose=true)
newton(F(10), D, 2, ftol=1e-8, verbose=true)

plot(x -> x + log(x), 0.1, 1; legend=false)
title!(L"f(x) = x + \ln(x)")
xlabel!(L"x")
ylabel!(L"f(x)")
xticks!(0.1:0.1:1)

function bisection(f, a, b; atol=1e-10, max_iter = Inf, verbose=false, exact_val=NaN)
    # Assume that a < b and f(a) * f(b) < 0 initially

    p = NaN
    iteration = 0

    if verbose
        @printf("%5s%15s", "k", "x_k")
        isnan(exact_val) || @printf("%30s", "Correct digits in frac. part")
        @printf("\n")
    end

    while iteration <= max_iter
        p_prev = p
        p = (a + b) / 2
        fp = f(p)
        fafp = f(a) * f(p)
        iteration += 1

        if verbose
            @printf("%5d%15f", iteration, p)
            print_correctdigits(p, exact_val)
        end

        (abs(p - p_prev) < atol || fafp == 0) && break

        fafp < 0 ? b = p : a = p
    end

    @printf("\n")
    return p
end


function fixedpoint(f, x_0; atol=1e-10, max_iter = Inf, verbose=false, exact_val=NaN)
    x_n = x_0
    iteration = 0

    if verbose
        @printf("%5s%15s", "k", "x_k")
        isnan(exact_val) || @printf("%30s", "Correct digits in frac. part")
        @printf("\n%5d%15f", iteration, x_n)
        print_correctdigits(x_n, exact_val)
    end

    while iteration <= max_iter
        x_prev = x_n
        x_n = f(x_n)
        iteration += 1

        if verbose
            @printf("%5d%15f", iteration, x_n)
            print_correctdigits(x_n, exact_val)
        end

        (abs(x_n - x_prev) < atol) && break
    end

    @printf("\n")
    return x_n
end


function secant(f, x_0, x_1; atol=1e-10, max_iter = Inf, verbose=false, exact_val=NaN)
    x_prev, x_n = x_0, x_1
    iteration = 1

    if verbose
        @printf("%5s%15s", "k", "x_k")
        isnan(exact_val) || @printf("%30s", "Correct digits in frac. part")
        @printf("\n%5d%15f", 0, x_0)
        print_correctdigits(x_0, exact_val)
        @printf("%5d%15f", 1, x_1)
        print_correctdigits(x_1, exact_val)
    end

    while iteration <= max_iter
        x_prev2 = x_prev
        x_prev = x_n
        x_n = x_prev - (f(x_prev) * (x_prev - x_prev2)) / (f(x_prev) - f(x_prev2))
        iteration += 1

        if verbose
            @printf("%5d%15f", iteration, x_n)
            print_correctdigits(x_n, exact_val)
        end

        (abs(x_n - x_prev) < atol) && break
    end

    @printf("\n")
    return x_n
end

    
# Root of x + ln(x), computed by WolframAlpha (for comparison purposes):
EXACT_ROOT = 0.5671432904097838729999686622103555497538

f = x -> x + log(x)
df = x -> 1 + 1/x
g = x-> exp(-x)

@printf("Bisection method result: x = %.15f\n\n\n", 
        bisection(f, 0.5, 0.6, verbose=true, exact_val=EXACT_ROOT))
@printf("Fixed point method result: x = %.15f\n\n\n", 
        fixedpoint(g, 0.5, verbose=true, exact_val=EXACT_ROOT))
@printf("Newton's method result: x = %.15f\n\n\n", 
        newton(f, df, 0.5, verbose=true, exact_val=EXACT_ROOT))
@printf("Secant method result: x = %.15f\n\n\n", 
        secant(f, 0.5, 0.6, verbose=true, exact_val=EXACT_ROOT))
