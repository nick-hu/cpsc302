
#=
CPSC 302 Assignment 8
Nicholas Hu
=#

using Plots
pyplot()
using LaTeXStrings
using LightGraphs, GraphLayout

u = 1:32
v = [1:30; 30; 32]
M = randn(32, 32)
Q, R = qr(M)
A = Q * diagm(u) * Q'
B = Q * diagm(v) * Q'

rayleigh = (A, v) -> (v' * A * v)[1]

eA, eB = maximum(eigvals(A)), maximum(eigvals(B))
v0 = randn(size(A, 1))
tol = 1e-10
err_A, err_B = [Inf], [Inf]

v = v0
while err_A[end] >= tol
    v = A*v
    v /= norm(v)
    push!(err_A, abs(rayleigh(A, v) - eA))
end

v = v0
while err_B[end] >= tol
    v = B*v
    v /= norm(v)
    push!(err_B, abs(rayleigh(B, v) - eB))
end

p = plot(err_A[2:end], yscale=:log10, label="Matrix A",
         title="Absolute error in power method vs. iteration number",
         xlabel="Iteration number",
         ylabel="Absolute error")
plot!(err_B[2:end], label="Matrix B")
display(p)
@printf("Ratio of iterations taken: %f\n", (length(err_A - 1)) / (length(err_B) - 1))

s = [1, 1, 2, 2, 3, 3, 3, 4, 5]
d = [2, 5, 3, 4, 4, 5, 6, 1, 1]
n = 6

G = DiGraph(n)
for e = 1:length(s)
    add_edge!(G, s[e], d[e])
end

A = full(adjacency_matrix(G))
loc_x, loc_y = layout_spring_adj(A)
draw_layout_adj(A, loc_x, loc_y, labels=collect(1:6), nodefillc="#FFFFFF")

H = Matrix{Float64}(n, n)

for i = 1:n
    d = outdegree(G, i)
    if d == 0
        H[i, :] = ones(n) / n
    else
        H[i, out_neighbors(G, i)] = 1 / d
    end
end

a = 0.85  # Damping factor
H2 = a * H + (1-a) * ones(n, n) / n
v = abs(real(eigvecs(H2')[:, 1]))  # Dominant eigenvector
v /= norm(v, 1)

bar(v, label="", xticks=1:n, ylims=[0, 1], title="PageRank of pages in web", 
    xlabel="Page number", ylabel="PageRank")
