using Distributions, StatsFuns

gicc(θ, ϕ, α, β) = logistic(α*ϕ / √(α^2+ϕ^2) * (θ - β))

using Plots
Plots.plot(-4:0.1:4, [gicc.(-4:0.1:4, ϕ, 2, 0) for ϕ in 0:0.1:4], line_z = (0:40)', c = :blues, colorbar = false)

