
#= 
CPSC 302 Assignment 1
Nicholas Hu
=#

# Centred difference approximation for derivative of f(x) = sin(x) at x0 = 1.2

using PyPlot

f(x) = sin(x)
df(x) = cos(x)
max_f3 = 1  # The maximum absolute value of the third derivative of f (-cos(x))

x0 = 1.2

h = logspace(0, -20, 41)
abs_err = abs(df(x0) - (f(x0 + h) - f(x0 - h)) ./ (2h))

@printf("%10s %20s\n", "h", "Absolute error")
for i in 1:length(h)
    @printf("%10.2e %20e\n", h[i], abs_err[i])
end

p1 = loglog(h, abs_err, marker="o", label="Absolute error")
p2 = loglog(h[1:12], (h[1:12].^2 / 6) * max_f3, linestyle="--", 
            label="Discretization error")  # Truncated plot for aesthetics
xlabel("h")
legend(loc="lower left")
title("Absolute and discretization errors in\n centred difference approximation")

u0 = log(3)
u_prev = u0

@printf("%5s %20s %20s\n", "n", "u_n", "u_n - u_(n-1)")
@printf("%5d %20f\n", 0, u0)

for n in 1:50
    u = 1/n - 0.5u_prev
    @printf("%5d %20f %20f\n", n, u, u - u_prev)
    u_prev = u
end
