clear; clc;

% ============================================================
% Compare periodic SCF solvers on the same assembled problem
% ============================================================
%
% This script:
%   1) Builds the same compact periodic test crystal
%   2) Runs the periodic triclinic Ewald polarization calculation with:
%         - direct
%         - matrix_iterative
%         - gs
%         - sor
%      for both:
%         - bare periodic operator
%         - Thole-corrected real-space operator
%   3) Compares all iterative methods against direct
%   4) Reports timing, convergence, residuals, energy errors, dipole errors
%
% Recommended use:
%   - validate GS/SOR implementation against direct
%   - compare efficiency of matrix_iterative vs GS/SOR
%   - tune omega for GS/SOR
%
% Notes:
%   - omega = 1.0 means pure GS
%   - 0 < omega < 1 means under-relaxed GS
%   - 1 < omega < 2 means true SOR
%
% Assumes:
%   - calc.run_polarization_calc supports params.scf.solver =:
%       'direct', 'matrix_iterative', 'gs', 'sor'
%   - GS/SOR implementation uses the assembled matrix T
%   - result struct contains:
%       res.mu
%       res.energy.total
%       res.T
%       res.Eext
%       res.sys
%       res.scf   (or compatible SCF metadata)
%

fprintf('\n============================================================\n');
fprintf('BUILD TEST SYSTEM\n');
fprintf('============================================================\n');

% -----------------------------
% Define compact crystal
% -----------------------------
crystal = struct();
crystal.cellpar = [6.5, 7.0, 6.2, 88.0, 101.0, 84.0];
crystal.lattice = [];

crystal.frac_coords = [
    0.10 0.15 0.20   % site 1: C1  (mol 1)
    0.16 0.18 0.24   % site 2: H1  (mol 1)
    0.42 0.20 0.28   % site 3: N1  (mol 2)
    0.48 0.24 0.31   % site 4: H2  (mol 2)
];

crystal.cart_coords = [];
crystal.mol_id = [1; 1; 2; 2];
crystal.site_label = {'C1'; 'H1'; 'N1'; 'H2'};
crystal.site_type  = {'C';  'H';  'N';  'H' };

% -----------------------------
% Define model
% -----------------------------
model = struct();
model.polarizable_types = {'C', 'N'};
model.alpha_by_type = struct();
model.alpha_by_type.C = 1.334;
model.alpha_by_type.N = 1.075;
model.thole_a = 0.39;

model.charge_pattern = struct();
model.charge_pattern.site_label = {'C1', 'H1'};
model.charge_pattern.delta_q = [0.10, -0.10];

% -----------------------------
% Neutral build for molecule lookup
% -----------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = true;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
res0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = res0.sys;

disp('Available molecules:')
disp(builder.list_molecules(sys0))

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

% ============================================================
% Solver settings to compare
% ============================================================
solvers = {
    struct('name','direct',           'solver','direct',           'omega', NaN)
    struct('name','matrix_iterative', 'solver','matrix_iterative', 'omega', 1.0)
    struct('name','gs',               'solver','gs',               'omega', 1.0)
    struct('name','sor_0p90',         'solver','sor',              'omega', 0.9)
    struct('name','sor_1p05',         'solver','sor',              'omega', 1.05)
};

tol = 1e-10;
maxIter = 5000;

% ============================================================
% Run two physics cases
% ============================================================
cases = {
    struct('name','bare',  'use_thole_real_space',false, 'thole_a',model.thole_a)
    struct('name','thole', 'use_thole_real_space',true,  'thole_a',model.thole_a)
};

allResults = struct();

for c = 1:numel(cases)
    caseInfo = cases{c};

    fprintf('\n============================================================\n');
    fprintf('CASE: %s\n', upper(caseInfo.name));
    fprintf('============================================================\n');

    for s = 1:numel(solvers)
        sol = solvers{s};

        params = util.default_params();
        params.ewald.mode = 'periodic_triclinic';
        params.ewald.auto = true;
        params.ewald.tol = 1e-8;
        params.ewald.boundary = 'tinfoil';
        params.ewald.rcut_fraction = 0.9;
        params.ewald.use_thole_real_space = caseInfo.use_thole_real_space;
        params.ewald.thole_a = caseInfo.thole_a;

        params.scf.solver = sol.solver;
        params.scf.tol = tol;
        params.scf.maxIter = maxIter;

        if ~isnan(sol.omega)
            params.scf.omega = sol.omega;
        end

        fprintf('Running %-16s ... ', sol.name);
        tStart = tic;
        res = calc.run_polarization_calc(crystal, model, opts, params);
        elapsed = toc(tStart);
        fprintf('done (%.4f s)\n', elapsed);

        out = struct();
        out.res = res;
        out.walltime = elapsed;
        out.solver_name = sol.name;
        out.solver = sol.solver;
        out.omega = sol.omega;

        % Robustly extract SCF info if present
        if isfield(res, 'scf') && ~isempty(res.scf)
            out.scf = res.scf;
        else
            out.scf = struct();
        end

        % Compute true linear-system residual from assembled T and returned mu
        out.relres_true = compute_true_relres_from_result(res);

        allResults.(caseInfo.name).(sol.name) = out;
    end
end

% ============================================================
% Report Ewald settings (from direct bare run)
% ============================================================
res_ref = allResults.bare.direct.res;

fprintf('\n============================================================\n');
fprintf('EWALD PARAMETERS USED\n');
fprintf('============================================================\n');
fprintf('alpha = %.16e\n', res_ref.Tmeta.alpha);
fprintf('rcut  = %.16e\n', res_ref.Tmeta.rcut);
fprintf('kcut  = %.16e\n', res_ref.Tmeta.kcut);
fprintf('num_kvec = %d\n', res_ref.Tmeta.num_kvec);
fprintf('real image bounds = [%d %d %d]\n', res_ref.Tmeta.real_image_bounds);

% ============================================================
% Compare against direct
% ============================================================
for c = 1:numel(cases)
    caseInfo = cases{c};
    ref = allResults.(caseInfo.name).direct.res;

    fprintf('\n============================================================\n');
    fprintf('COMPARISON AGAINST DIRECT: %s\n', upper(caseInfo.name));
    fprintf('============================================================\n');

    fprintf('%-18s %-10s %-8s %-12s %-14s %-14s %-14s %-14s %-14s\n', ...
        'solver', 'conv', 'iters', 'time(s)', 'relres_scf', 'relres_true', ...
        '|dU|', 'max|dmu_i|', '||dmu||_F');
    fprintf('%s\n', repmat('-',1,124));

    for s = 1:numel(solvers)
        sol = solvers{s};
        out = allResults.(caseInfo.name).(sol.name);
        res = out.res;

        [convFlag, nIter] = unpack_scf_info(out.scf, sol.solver);

        solverRelres = NaN;
        if isfield(out.scf, 'relres') && ~isempty(out.scf.relres)
            solverRelres = out.scf.relres;
        end

        dU = abs(res.energy.total - ref.energy.total);
        dmu = res.mu - ref.mu;
        maxdmu = max(sqrt(sum(dmu.^2,2)));
        frodmu = norm(dmu, 'fro');

        fprintf('%-18s %-10d %-8d %-12.4f %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', ...
            sol.name, convFlag, nIter, out.walltime, solverRelres, out.relres_true, ...
            dU, maxdmu, frodmu);
    end
end

% ============================================================
% Optional omega scan for SOR
% ============================================================
fprintf('\n============================================================\n');
fprintf('OPTIONAL OMEGA SCAN (THOLE CASE)\n');
fprintf('============================================================\n');

omegaList = [0.7 0.8 0.9 1.0 1.05 1.1 1.2];
omegaResults = struct([]);

ref = allResults.thole.direct.res;

fprintf('%-10s %-10s %-8s %-12s %-14s %-14s %-14s\n', ...
    'omega', 'conv', 'iters', 'time(s)', 'relres_true', '|dU|', 'max|dmu_i|');
fprintf('%s\n', repmat('-',1,90));

for k = 1:numel(omegaList)
    omega = omegaList(k);

    params = util.default_params();
    params.ewald.mode = 'periodic_triclinic';
    params.ewald.auto = true;
    params.ewald.tol = 1e-8;
    params.ewald.boundary = 'tinfoil';
    params.ewald.rcut_fraction = 0.9;
    params.ewald.use_thole_real_space = true;
    params.ewald.thole_a = model.thole_a;

    params.scf.solver = 'sor';
    params.scf.omega = omega;
    params.scf.tol = tol;
    params.scf.maxIter = maxIter;

    tStart = tic;
    res = calc.run_polarization_calc(crystal, model, opts, params);
    elapsed = toc(tStart);

    relres_true = compute_true_relres_from_result(res);
    [convFlag, nIter] = unpack_scf_info(get_scf_struct(res), 'sor');

    dU = abs(res.energy.total - ref.energy.total);
    dmu = res.mu - ref.mu;
    maxdmu = max(sqrt(sum(dmu.^2,2)));

    fprintf('%-10.3f %-10d %-8d %-12.4f %-14.6e %-14.6e %-14.6e\n', ...
        omega, convFlag, nIter, elapsed, relres_true, dU, maxdmu);

    omegaResults(k).omega = omega; %#ok<SAGROW>
    omegaResults(k).res = res; %#ok<SAGROW>
    omegaResults(k).time = elapsed; %#ok<SAGROW>
    omegaResults(k).relres_true = relres_true; %#ok<SAGROW>
    omegaResults(k).dU = dU; %#ok<SAGROW>
    omegaResults(k).maxdmu = maxdmu; %#ok<SAGROW>
end

fprintf('\nDone.\n');

% ============================================================
% Local helper functions
% ============================================================

function relres = compute_true_relres_from_result(res)
% Compute relative residual for:
%   (I - A*T) mu = A*Eext
%
% using returned res.sys, res.T, res.Eext, res.mu

    sys = res.sys;
    T = res.T;
    mu = res.mu;
    Eext = res.Eext;

    nSites = size(mu,1);
    polMask = logical(sys.site_is_polarizable(:));
    alpha = sys.site_alpha(:);

    muVec = reshape(mu.', [], 1);     % [mux1 muy1 muz1 mux2 ...]^T
    EVec  = reshape(Eext.', [], 1);

    Avec = zeros(3*nSites, 1);
    for i = 1:nSites
        Ii = (3*i-2):(3*i);
        if polMask(i)
            Avec(Ii) = alpha(i) * EVec(Ii);
        end
    end

    r = Avec - (muVec - apply_A_to_Tmu(T, muVec, alpha, polMask));
    denom = norm(Avec);
    if denom == 0
        denom = 1.0;
    end
    relres = norm(r) / denom;
end

function ATmu = apply_A_to_Tmu(T, muVec, alpha, polMask)
% Return A*(T*mu), where A is block diagonal with alpha_i*I3 on polarizable sites.

    y = T * muVec;
    nSites = numel(alpha);
    ATmu = zeros(size(y));

    for i = 1:nSites
        Ii = (3*i-2):(3*i);
        if polMask(i)
            ATmu(Ii) = alpha(i) * y(Ii);
        end
    end
end

function scf = get_scf_struct(res)
    if isfield(res, 'scf') && ~isempty(res.scf)
        scf = res.scf;
    else
        scf = struct();
    end
end

function [convFlag, nIter] = unpack_scf_info(scf, solverName)
    convFlag = NaN;
    nIter = NaN;

    if strcmpi(solverName, 'direct')
        convFlag = 1;
        nIter = 1;
        return;
    end

    if isfield(scf, 'converged') && ~isempty(scf.converged)
        convFlag = scf.converged;
    end

    if isfield(scf, 'nIter') && ~isempty(scf.nIter)
        nIter = scf.nIter;
    elseif isfield(scf, 'iters') && ~isempty(scf.iters)
        nIter = scf.iters;
    end
end