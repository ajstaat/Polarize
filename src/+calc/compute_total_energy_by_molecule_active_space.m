function tbl = compute_total_energy_by_molecule_active_space(sys, problem, mu, Eext, Tpol)
%COMPUTE_TOTAL_ENERGY_BY_MOLECULE_ACTIVE_SPACE Molecule-resolved energy table.

    mu_pol = mu(problem.activeSites, :);
    Eext_pol = Eext(problem.activeSites, :);
    mol_pol = sys.site_mol_id(problem.activeSites);

    mu_pol_vec = util.stack_xyz(mu_pol);
    Eext_pol_vec = util.stack_xyz(Eext_pol);
    Tmu_pol_vec = Tpol * mu_pol_vec;

    uniqueMol = unique(mol_pol);
    nMol = numel(uniqueMol);

    ext_energy = zeros(nMol,1);
    dd_energy  = zeros(nMol,1);

    for m = 1:nMol
        mol = uniqueMol(m);
        idx_sites = find(mol_pol == mol);

        vec_idx = zeros(3*numel(idx_sites),1);
        for k = 1:numel(idx_sites)
            a = idx_sites(k);
            vec_idx(3*k-2:3*k) = (3*a-2):(3*a);
        end

        ext_energy(m) = -0.5 * dot(mu_pol_vec(vec_idx), Eext_pol_vec(vec_idx));
        dd_energy(m)  = -0.5 * dot(mu_pol_vec(vec_idx), Tmu_pol_vec(vec_idx));
    end

    tbl = table(uniqueMol, ext_energy, dd_energy, ext_energy + dd_energy, ...
        'VariableNames', {'mol_id', 'external_charge_dipole', 'dipole_dipole', 'total'});
end