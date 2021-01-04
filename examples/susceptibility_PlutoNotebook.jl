### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 66923a50-0fa8-11eb-1547-cf8e289818b0
using Revise

# ╔═╡ 8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
using Pkg

# ╔═╡ 93125894-0fa8-11eb-0ff6-211feff97889
using PolaronMobility

# ╔═╡ f0529460-0fa8-11eb-1bec-1feb30645d42
using Plots

# ╔═╡ 8ae04a5a-0fa8-11eb-2821-8dc46ee73269
Pkg.develop("PolaronMobility")

# ╔═╡ 386fd0a4-0faa-11eb-14e3-a3cdf3f3cb5d
α=5.0

# ╔═╡ a86c04ae-0faa-11eb-3565-7bd3108906b2
feynmanvw(α) # Variational solution, athermal (original) action

# ╔═╡ 3e91bf88-0faa-11eb-33db-8b7fd24cf1eb
βred=5

# ╔═╡ b3e10abe-0faa-11eb-308e-91c4c0886b8b
feynmanvw(α, βred) # Variational solution, thermal (Osaka) action

# ╔═╡ 4306f2f4-0faa-11eb-1d32-dd8ade1d930a
ω=2.25E12

# ╔═╡ 4b26ec8c-0faa-11eb-2af1-e3f5e39fb358
meff=0.12

# ╔═╡ 34a7c84e-0faa-11eb-08b6-83f4ab86b689
nurange=0:0.1:22

# ╔═╡ ad10a50c-0fa8-11eb-30a9-1b1bb4fc26b7
X=PolaronMobility.ImX(nurange,feynmanvw(α, βred)...,βred,α,ω,meff*PolaronMobility.MassElectron)

# ╔═╡ f384207c-0fa8-11eb-29b6-9b35b33aa526
plot(X.nu, X.μ, label="Mobility")

# ╔═╡ a3cb113e-0fa9-11eb-1093-d7f12e1cf14b
plot(X.nu, X.ImX, label="ImX")

# ╔═╡ Cell order:
# ╠═66923a50-0fa8-11eb-1547-cf8e289818b0
# ╠═8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
# ╠═8ae04a5a-0fa8-11eb-2821-8dc46ee73269
# ╠═93125894-0fa8-11eb-0ff6-211feff97889
# ╠═386fd0a4-0faa-11eb-14e3-a3cdf3f3cb5d
# ╠═a86c04ae-0faa-11eb-3565-7bd3108906b2
# ╠═3e91bf88-0faa-11eb-33db-8b7fd24cf1eb
# ╠═b3e10abe-0faa-11eb-308e-91c4c0886b8b
# ╠═4306f2f4-0faa-11eb-1d32-dd8ade1d930a
# ╠═4b26ec8c-0faa-11eb-2af1-e3f5e39fb358
# ╠═34a7c84e-0faa-11eb-08b6-83f4ab86b689
# ╠═ad10a50c-0fa8-11eb-30a9-1b1bb4fc26b7
# ╠═f0529460-0fa8-11eb-1bec-1feb30645d42
# ╠═f384207c-0fa8-11eb-29b6-9b35b33aa526
# ╠═a3cb113e-0fa9-11eb-1093-d7f12e1cf14b
