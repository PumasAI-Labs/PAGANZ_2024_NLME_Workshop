using AlgebraOfGraphics, CairoMakie, ForwardDiff, Optim, LaTeXStrings

f(x) = (-x * (x - 2) * (x - 4) * (x - 7.5) + 30) / 200
x = 0:0.01:8.1
mode = Optim.optimize(x -> -f(x[1]), [6.0]).minimizer[1]
h = ForwardDiff.hessian(x -> f(x[1]), [mode])[1]
g(x) = f(mode) + h / 2 * (x - mode)^2
firstind = findfirst(x -> g(x) > -0.5, x)
gx = x[firstind:end]

fig = Figure()
for (col, j) in enumerate(1:4:16)
    ax = Axis(fig[1,col], title = latexstring("\$$(j == 1 ? "" : "$j \\cdot ")\\log f(\\eta)\$"))
    line1 = CairoMakie.lines!(ax, x, j .* f.(x))
    line2 = CairoMakie.lines!(ax, gx, j .* g.(gx), color = :red, label = "Laplace approximation")
    vlines!(ax, [mode], color = :black, linestyle = :dash)
    ax = Axis(fig[2,col], title = latexstring("\$f(\\eta)$(j == 1 ? "" : "^{$j}")\$"))
    line1 = CairoMakie.lines!(ax, x, exp.(j .* f.(x)))
    line2 = CairoMakie.lines!(ax, gx, exp.(j .* g.(gx)), color = :red, label = "Laplace approximation")
    vlines!(ax, [mode], color = :black, linestyle = :dash)
end
fig

save("supporting_material/laplace_method.png", fig)
