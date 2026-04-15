function out = find_enclosing_supercell_for_molecules(filename, varargin)
%FIND_ENCLOSING_SUPERCELL_FOR_MOLECULES
% For each unwrapped unit-cell molecule from a CONTCAR/POSCAR, find a
% lattice-aligned supercell box that fully encloses that molecule as a
% contiguous object.
%
% This does NOT build a Polarize system. It only analyzes whole molecules
% reconstructed from the source file.
%
% Output struct fields:
%   .molecules        cell array of unwrapped molecules from io.unwrap_all_contcar_molecules
%   .table            table with one row per base molecule:
%                       base_mol_id
%                       n_sites
%                       frac_min
%                       frac_max
%                       frac_shift
%                       dims
%                       volume_cells
%   .best             struct for the smallest-volume choice
%
% Example:
%   B = builder.find_enclosing_supercell_for_molecules(filename, ...
%       'BondScale', 1.20, 'SortMolecules', false, 'Buffer', 0.10);

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Buffer', 0.10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, filename, varargin{:});
    opt = p.Results;

    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    nMol = numel(mols);

    base_mol_id = zeros(nMol,1);
    n_sites = zeros(nMol,1);
    frac_min = zeros(nMol,3);
    frac_max = zeros(nMol,3);
    frac_shift = zeros(nMol,3);
    dims = zeros(nMol,3);
    volume_cells = zeros(nMol,1);

    for m = 1:nMol
        Mol = mols{m};
        f = Mol.fracUnwrapped;

        % Shift so the molecule sits just inside the box with a small buffer
        fmin0 = min(f, [], 1);
        shift = opt.Buffer - fmin0;
        fshift = f + shift;

        fmin = min(fshift, [], 1);
        fmax = max(fshift, [], 1);

        boxDims = ceil(fmax);

        base_mol_id(m) = m;
        n_sites(m) = size(f,1);
        frac_min(m,:) = fmin;
        frac_max(m,:) = fmax;
        frac_shift(m,:) = shift;
        dims(m,:) = boxDims;
        volume_cells(m) = prod(boxDims);
    end

    T = table(base_mol_id, n_sites, frac_min, frac_max, frac_shift, dims, volume_cells);

    [~, ibest] = min(volume_cells);

    best = struct();
    best.base_mol_id = T.base_mol_id(ibest);
    best.n_sites = T.n_sites(ibest);
    best.frac_min = T.frac_min(ibest,:);
    best.frac_max = T.frac_max(ibest,:);
    best.frac_shift = T.frac_shift(ibest,:);
    best.dims = T.dims(ibest,:);
    best.volume_cells = T.volume_cells(ibest);
    best.molecule = mols{best.base_mol_id};

    out = struct();
    out.molecules = mols;
    out.table = T;
    out.best = best;
end