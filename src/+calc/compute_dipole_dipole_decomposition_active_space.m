function decomp = compute_dipole_dipole_decomposition_active_space(sys, problem, mu, Tpol)
%COMPUTE_DIPOLE_DIPOLE_DECOMPOSITION_ACTIVE_SPACE Pairwise dd decomposition.

    nPol = problem.nPolSites;
    mu_pol = mu(problem.activeSites, :);

    pair_i = [];
    pair_j = [];
    pair_energy = [];

    for a = 1:nPol
        Ia = (3*a-2):(3*a);
        mua = mu_pol(a,:).';
        for b = (a+1):nPol
            Ib = (3*b-2):(3*b);
            mub = mu_pol(b,:).';

            eab = -mua.' * Tpol(Ia, Ib) * mub;

            pair_i(end+1,1) = problem.activeSites(a); %#ok<AGROW>
            pair_j(end+1,1) = problem.activeSites(b); %#ok<AGROW>
            pair_energy(end+1,1) = eab; %#ok<AGROW>
        end
    end

    decomp = struct();
    decomp.site_i = pair_i;
    decomp.site_j = pair_j;
    decomp.energy = pair_energy;
    decomp.total = sum(pair_energy);
end