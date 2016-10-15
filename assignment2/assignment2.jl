
#= 
CPSC 302 Assignment 2
Nicholas Hu
=#

using PyPlot

function horner_eval(coeffs, x)
    y = 0
    
    for c in coeffs
        y = y .* x + c
    end
    
    return y
end

f(x) = (x-2).^9
f_horner(x) = horner_eval([1, -18, 144, -672, 2016, -4032, 5376, -4608, 2304, -512], x)

x = 1.94:0.001:2.08

plot(x, f(x), "r")
axis("tight")
title(L"$f(x) = (x-2)^9$ evaluated directly")
xlabel(L"$x$")
ylabel(L"$f(x)$")

figure()

plot(x, f_horner(x), "b")
axis("tight")
title(L"$f(x) = (x-2)^9$ evaluated using Horner's rule")
xlabel(L"$x$")
ylabel(L"$f(x)$")
