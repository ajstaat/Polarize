%% find_quasi_isotropic_supercells
%
% Geometry-only helper for choosing periodic finite-size convergence cells.
%
% Goal:
%   Enumerate integer supercell replications [na nb nc] and rank them by
%   quasi-isotropic growth around the charged pair / periodic images.
%
% Key metrics:
%
%   Hsuper rows are direct lattice vectors:
%
%       cart = frac * H
%
%   Lmin:
%       shortest nonzero lattice translation length.
%
%   R_eff:
%       0.5 * Lmin, the radius of the largest sphere guaranteed not to
%       overlap its nearest periodic image under the minimum-image picture.
%
%   aspect_lengths:
%       max(row-vector length) / min(row-vector length)
%
%   aspect_image:
%       rough anisotropy of the shortest-image shell.
%
%   p3m_mesh:
%       suggested mesh preserving approximately the target real-space mesh
%       spacing.
%
% This script does not build the full molecular system and does not run
% polarization. It is intended to choose a clean convergence ladder before
% launching expensive workflows.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Find quasi-isotropic supercells for finite-size convergence\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% User controls
% -------------------------------------------------------------------------

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Search ranges.
%
% Adjust these upward if you want to consider larger cells.
cfg.na_range = 2:8;
cfg.nb_range = 4:20;
cfg.nc_range = 2:8;

% Primitive template info.
cfg.n_sites_primitive = 288;
cfg.n_molecules_primitive = 4;

% Practical caps. These are geometry-screening caps, not hard code limits.
cfg.max_sites = 150000;
cfg.max_volume_factor = inf;

% Current known/reference mesh choices.
%
% If you consider [3 11 3] with [64 216 64] spacing-consistent, this sets
% the target grid spacing. The workflow will suggest meshes for every
% candidate using approximately the same spacing.
cfg.reference_supercell = [3 11 3];
cfg.reference_mesh = [64 216 64];

% Optional FFT-friendly rounding.
%
% 'none'      : use raw rounded mesh.
% 'even'      : round up to next even integer.
% 'smooth235' : round up to a number whose prime factors are only 2,3,5.
cfg.mesh_rounding = 'smooth235';

% Minimum mesh count per direction.
cfg.min_mesh = [16 16 16];

% Ranking controls.
%
% Score is lower-is-better:
%   score = wAspect * log(aspect_lengths)^2
%         + wImage  * log(aspect_image)^2
%         + wSize   * volume_penalty
%
% For choosing a convergence ladder, we mostly want large R_eff with decent
% shape, not necessarily the smallest score globally.
cfg.wAspect = 1.0;
cfg.wImage = 0.5;
cfg.wSize = 0.02;

% Recommended ladder settings.
cfg.n_ladder = 10;
cfg.min_R_spacing_bohr = 5.0;

% Print controls.
cfg.print_top_n = 40;
cfg.print_ladder_n = 12;

%% ------------------------------------------------------------------------
% Read base lattice
% -------------------------------------------------------------------------

if ~isfile(cfg.filename)
    error('Input file not found:\n  %s', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

H0 = local_get_lattice_rows(crystal);

fprintf('\nPrimitive lattice H rows / bohr:\n');
disp(H0);

baseLengths = vecnorm(H0, 2, 2);
baseVolume = abs(det(H0));
baseLmin = local_shortest_lattice_translation_bruteforce(H0, 3);

fprintf('Primitive row-vector lengths / bohr: [%.6f %.6f %.6f]\n', ...
    baseLengths(1), baseLengths(2), baseLengths(3));
fprintf('Primitive volume / bohr^3        : %.8e\n', baseVolume);
fprintf('Primitive shortest translation   : %.8f bohr\n', baseLmin);

%% ------------------------------------------------------------------------
% Reference mesh spacing
% -------------------------------------------------------------------------

Href = diag(cfg.reference_supercell) * H0;
refLengths = vecnorm(Href, 2, 2);
targetSpacing = refLengths ./ cfg.reference_mesh(:);

fprintf('\nReference spacing from supercell [%d %d %d], mesh [%d %d %d]:\n', ...
    cfg.reference_supercell(1), cfg.reference_supercell(2), cfg.reference_supercell(3), ...
    cfg.reference_mesh(1), cfg.reference_mesh(2), cfg.reference_mesh(3));
fprintf('  target spacing / bohr = [%.6f %.6f %.6f]\n', ...
    targetSpacing(1), targetSpacing(2), targetSpacing(3));

%% ------------------------------------------------------------------------
% Enumerate candidates
% -------------------------------------------------------------------------

rows = [];

for na = cfg.na_range
    for nb = cfg.nb_range
        for nc = cfg.nc_range
            rep = [na nb nc];

            nSites = cfg.n_sites_primitive * prod(rep);
            nMols = cfg.n_molecules_primitive * prod(rep);
            volFactor = prod(rep);

            if nSites > cfg.max_sites
                continue;
            end

            if volFactor > cfg.max_volume_factor
                continue;
            end

            H = diag(rep) * H0;

            lengths = vecnorm(H, 2, 2);
            volume = abs(det(H));

            Lmin = local_shortest_lattice_translation_bruteforce(H, 3);
            R_eff = 0.5 * Lmin;

            aspectLengths = max(lengths) / min(lengths);

            shell = local_shortest_shell_metrics(H, 3);
            aspectImage = shell.max_shell_norm / max(shell.min_shell_norm, eps);

            rawMesh = ceil(lengths ./ targetSpacing);
            rawMesh = max(rawMesh(:).', cfg.min_mesh);

            mesh = local_round_mesh(rawMesh, cfg.mesh_rounding);

            meshSpacing = lengths(:).' ./ mesh;
            meshSpacingRatio = max(meshSpacing ./ targetSpacing(:).') / ...
                min(meshSpacing ./ targetSpacing(:).');

            volumePenalty = log(max(volFactor, 1))^2;

            score = cfg.wAspect * log(aspectLengths)^2 + ...
                cfg.wImage * log(aspectImage)^2 + ...
                cfg.wSize * volumePenalty;

            row = struct();
            row.na = na;
            row.nb = nb;
            row.nc = nc;
            row.nSites = nSites;
            row.nMolecules = nMols;
            row.volumeFactor = volFactor;
            row.volume_bohr3 = volume;

            row.La = lengths(1);
            row.Lb = lengths(2);
            row.Lc = lengths(3);

            row.Lmin = Lmin;
            row.R_eff = R_eff;
            row.aspect_lengths = aspectLengths;
            row.aspect_image = aspectImage;
            row.n_shortest_images = shell.n_shortest;

            row.mesh_a = mesh(1);
            row.mesh_b = mesh(2);
            row.mesh_c = mesh(3);
            row.mesh_points = prod(mesh);
            row.mesh_spacing_a = meshSpacing(1);
            row.mesh_spacing_b = meshSpacing(2);
            row.mesh_spacing_c = meshSpacing(3);
            row.mesh_spacing_ratio = meshSpacingRatio;

            row.score = score;

            rows = [rows; row]; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    error('No candidate supercells found. Relax search ranges or max_sites.');
end

T = struct2table(rows);

% Useful derived labels.
T.rep = strings(height(T), 1);
T.mesh = strings(height(T), 1);

for i = 1:height(T)
    T.rep(i) = sprintf('[%d %d %d]', T.na(i), T.nb(i), T.nc(i));
    T.mesh(i) = sprintf('[%d %d %d]', T.mesh_a(i), T.mesh_b(i), T.mesh_c(i));
end

%% ------------------------------------------------------------------------
% Print top candidates by score
% -------------------------------------------------------------------------

Tscore = sortrows(T, {'score', 'R_eff'}, {'ascend', 'descend'});

fprintf('\n============================================================\n');
fprintf('Top quasi-isotropic candidates by score\n');
fprintf('============================================================\n');

local_print_candidate_table(Tscore, min(cfg.print_top_n, height(Tscore)));

%% ------------------------------------------------------------------------
% Print candidates by R_eff, with reasonable shape filters
% -------------------------------------------------------------------------

shapeMask = T.aspect_lengths <= 2.5 & T.aspect_image <= 1.5;
Tradius = sortrows(T(shapeMask, :), {'R_eff', 'score'}, {'ascend', 'ascend'});

fprintf('\n============================================================\n');
fprintf('Radius-ordered candidates with moderate shape filters\n');
fprintf('  filters: aspect_lengths <= 2.5, aspect_image <= 1.5\n');
fprintf('============================================================\n');

if isempty(Tradius)
    fprintf('No candidates passed moderate shape filters.\n');
else
    local_print_candidate_table(Tradius, min(cfg.print_top_n, height(Tradius)));
end

%% ------------------------------------------------------------------------
% Build recommended convergence ladder
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Recommended convergence ladder\n');
fprintf('============================================================\n');

if isempty(Tradius)
    fprintf('Using all candidates because no candidate passed shape filters.\n');
    Tpool = sortrows(T, {'R_eff', 'score'}, {'ascend', 'ascend'});
else
    Tpool = Tradius;
end

ladderIdx = local_select_ladder(Tpool, cfg.n_ladder, cfg.min_R_spacing_bohr);
Tladder = Tpool(ladderIdx, :);

local_print_candidate_table(Tladder, min(cfg.print_ladder_n, height(Tladder)));

fprintf('\nSuggested workflow entries:\n\n');

for i = 1:height(Tladder)
    fprintf('%% R_eff = %.3f bohr, aspect = %.3f, nSites = %d\n', ...
        Tladder.R_eff(i), Tladder.aspect_lengths(i), Tladder.nSites(i));
    fprintf('cfg.supercellSize = [%d %d %d];\n', ...
        Tladder.na(i), Tladder.nb(i), Tladder.nc(i));
    fprintf('cfg.p3m.mesh_size = [%d %d %d];\n\n', ...
        Tladder.mesh_a(i), Tladder.mesh_b(i), Tladder.mesh_c(i));
end

%% ------------------------------------------------------------------------
% Compare current / known cells explicitly
% -------------------------------------------------------------------------

knownCells = [
    3  5 3
    3 11 3
    5 13 5
    6 16 4
];

fprintf('\n============================================================\n');
fprintf('Known cells from recent runs\n');
fprintf('============================================================\n');

knownMask = false(height(T), 1);

for k = 1:size(knownCells, 1)
    rep = knownCells(k, :);
    knownMask = knownMask | ...
        (T.na == rep(1) & T.nb == rep(2) & T.nc == rep(3));
end

Tknown = sortrows(T(knownMask, :), {'R_eff'}, {'ascend'});

if isempty(Tknown)
    fprintf('None of the known cells are in the search range.\n');
else
    local_print_candidate_table(Tknown, height(Tknown));
end

fprintf('\nDone.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function H = local_get_lattice_rows(crystal)
% Return direct lattice vectors as rows in bohr.

if isfield(crystal, 'lattice')
    H = crystal.lattice;
elseif isfield(crystal, 'H')
    H = crystal.H;
elseif isfield(crystal, 'cell') && isfield(crystal.cell, 'H')
    H = crystal.cell.H;
elseif isfield(crystal, 'lattice_vectors')
    H = crystal.lattice_vectors;
else
    error('Could not find lattice matrix in crystal struct.');
end

H = double(H);

if ~isequal(size(H), [3 3])
    error('Lattice matrix must be 3x3.');
end
end

function Lmin = local_shortest_lattice_translation_bruteforce(H, searchReach)
% Brute-force shortest nonzero lattice translation for row-vector lattice.
%
% Translation vector:
%   t = n * H
%
% where n is integer row vector.

if nargin < 2
    searchReach = 3;
end

Lmin = inf;

for i = -searchReach:searchReach
    for j = -searchReach:searchReach
        for k = -searchReach:searchReach
            n = [i j k];

            if all(n == 0)
                continue;
            end

            t = n * H;
            d = norm(t);

            if d < Lmin
                Lmin = d;
            end
        end
    end
end

if ~isfinite(Lmin)
    error('Failed to compute shortest lattice translation.');
end
end

function shell = local_shortest_shell_metrics(H, searchReach)
% Compute approximate shortest-shell anisotropy from brute-force image list.

if nargin < 2
    searchReach = 3;
end

norms = [];

for i = -searchReach:searchReach
    for j = -searchReach:searchReach
        for k = -searchReach:searchReach
            n = [i j k];

            if all(n == 0)
                continue;
            end

            t = n * H;
            norms(end+1, 1) = norm(t); %#ok<AGROW>
        end
    end
end

norms = sort(norms);
dmin = norms(1);

tol = max(1e-8, 1e-8 * dmin);
shortest = norms(abs(norms - dmin) <= tol);

% Also inspect the first few low-lying image lengths to get a rough sense of
% near-shell anisotropy. This is not a rigorous Wigner-Seitz shape metric,
% but it catches obviously skinny cells.
nInspect = min(6, numel(norms));
low = norms(1:nInspect);

shell = struct();
shell.min_shell_norm = dmin;
shell.max_shell_norm = max(low);
shell.n_shortest = numel(shortest);
shell.low_norms = low;
end

function mesh = local_round_mesh(rawMesh, mode)
rawMesh = ceil(double(rawMesh(:).'));

switch lower(char(string(mode)))
    case 'none'
        mesh = rawMesh;

    case 'even'
        mesh = 2 * ceil(rawMesh / 2);

    case 'smooth235'
        mesh = zeros(1, 3);

        for i = 1:3
            mesh(i) = local_next_smooth235(rawMesh(i));
        end

    otherwise
        error('Unknown mesh rounding mode "%s".', mode);
end
end

function n = local_next_smooth235(n0)
% Small helper for FFT-friendly-ish mesh sizes.
%
% Returns the smallest integer >= n0 whose prime factors are only 2, 3, 5.

n0 = max(1, ceil(n0));
n = n0;

while true
    m = n;

    for p = [2 3 5]
        while mod(m, p) == 0
            m = m / p;
        end
    end

    if m == 1
        return;
    end

    n = n + 1;
end
end

function idx = local_select_ladder(Tpool, nMax, minSpacing)
idx = [];

if isempty(Tpool)
    return;
end

lastR = -inf;

for i = 1:height(Tpool)
    R = Tpool.R_eff(i);

    if isempty(idx) || (R - lastR) >= minSpacing
        idx(end+1, 1) = i; %#ok<AGROW>
        lastR = R;
    end

    if numel(idx) >= nMax
        break;
    end
end

% If spacing was too strict and produced too few points, fill by largest
% remaining radius gaps / low score.
if numel(idx) < min(nMax, height(Tpool))
    missing = setdiff((1:height(Tpool)).', idx);
    Tmiss = Tpool(missing, :);
    [~, order] = sortrows(Tmiss, {'R_eff', 'score'}, {'ascend', 'ascend'});
    fill = missing(order);

    for k = 1:numel(fill)
        idx(end+1, 1) = fill(k); %#ok<AGROW>

        if numel(idx) >= min(nMax, height(Tpool))
            break;
        end
    end

    idx = sort(idx);
end
end

function local_print_candidate_table(T, n)
n = min(n, height(T));

if n == 0
    fprintf('No candidates.\n');
    return;
end

cols = {'rep', 'nSites', 'R_eff', 'Lmin', ...
    'La', 'Lb', 'Lc', ...
    'aspect_lengths', 'aspect_image', ...
    'mesh', 'mesh_points', 'mesh_spacing_ratio', 'score'};

disp(T(1:n, cols));
end