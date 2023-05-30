ρ_array = 0.1:0.1:30.0
E_s = []
gauss_para = GaussParameter(100)
for ρ in ρ_array
    ge = GreensElement(0.0, 0.0, 5.0, 6.0, ρ, 10.0, 1.0, 1e-4)
    Es = QuaisEwald_Es_pair(1.0, 1.0, 1.0, ge, gauss_para)
    push!(E_s, Es)
end
plot(log10.(ρ_array), log10.(abs.(E_s)))

∂E = [(log10(abs(E_s[i + 10])) - log10(abs(E_s[i])))/(log10(ρ_array[i + 10]) - log10(ρ_array[i])) for i in 1:length(ρ_array) - 10]
x = [log10(ρ_array[i + 5])  for i in 1:length(ρ_array) - 10]

plot(x, ∂E, ylim = [-3, 3])

