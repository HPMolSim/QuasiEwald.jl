# benchmark the summations in the long range interaction
using BenchmarkTools, Profile, StatProfilerHTML
ge = GreensElement( - 0.8, 0.9, 10.0, 1e-4)
kx = rand()
ky = rand()
k_set = (kx, ky, sqrt(kx^2 + ky^2))
n_atoms = 1000
q = rand(n_atoms);
coords = [Point(rand(), rand(), rand()) for i in 1:n_atoms];
z_coords = [coords[i][3] for i in 1:n_atoms];
z_list = sortperm(z_coords);
@benchmark energy_k_sum($k_set, $q, $coords, $z_list, $ge)
@profilehtml (for i in 1:100000 energy_k_sum(k_set, q, coords, z_list, ge) end)
