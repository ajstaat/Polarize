%% tutorial_02_vasp_io
% This tutorial demonstrates the first IO layer in Polarize:
%
%   1. reading a minimal VASP POSCAR/CONTCAR-style file
%   2. inspecting raw parser output
%   3. converting a system-like struct to internal atomic units
%
% This tutorial does not yet unwrap molecules or build bond graphs.

clear; clc;

fprintf('\n============================================================\n');
fprintf('Tutorial 02: VASP IO and unit normalization\n');
fprintf('============================================================\n');

%% 1. Create a tiny temporary POSCAR-style file

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() local_cleanup(tmpDir));

filename = fullfile(tmpDir, 'POSCAR_tiny');

fid = fopen(filename, 'w');
if fid < 0
    error('Failed to open temporary POSCAR file.');
end

fprintf(fid, 'Tiny tutorial POSCAR\n');
fprintf(fid, '1.0\n');
fprintf(fid, '10.0 0.0 0.0\n');
fprintf(fid, '0.0 11.0 0.0\n');
fprintf(fid, '0.0 0.0 12.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '1 2\n');
fprintf(fid, 'Direct\n');
fprintf(fid, '0.0 0.0 0.0\n');
fprintf(fid, '0.5 0.0 0.0\n');
fprintf(fid, '0.0 0.5 0.0\n');

fclose(fid);

fprintf('\nWrote temporary POSCAR:\n  %s\n', filename);

%% 2. Parse the file

S = io.read_vasp_structure(filename);

fprintf('\nRaw VASP parser output:\n');
fprintf('  comment    = %s\n', S.comment);
fprintf('  scale      = %.6f\n', S.scale);
fprintf('  coord_type = %s\n', S.coord_type);
fprintf('  natoms     = %d\n', S.natoms);
fprintf('  lattice    = [%d x %d]\n', size(S.lattice,1), size(S.lattice,2));
fprintf('  frac       = [%d x %d]\n', size(S.frac,1), size(S.frac,2));
fprintf('  cart       = [%d x %d]\n', size(S.cart,1), size(S.cart,2));

fprintf('\nSpecies labels:\n');
disp(S.species(:).');

fprintf('\nCoordinate convention:\n');
fprintf('  For Direct VASP coordinates, parser computes cart = frac * lattice.\n');
fprintf('  frac/cart roundtrip error = %.3e\n', ...
    norm(S.frac - S.cart / S.lattice, 'fro'));

%% 3. Convert a system-like struct to atomic units

sys = struct();
sys.site_pos = S.cart;             % currently Angstrom-like for this tutorial
sys.site_alpha = [1.750; 0.696; 0.696];  % Angstrom^3 example values
sys.site_charge = zeros(S.natoms, 1);

sys.units.length = 'angstrom';
sys.units.alpha = 'angstrom^3';
sys.units.charge = 'elementary_charge';

sysAU = io.convert_to_atomic_units(sys);
io.assert_atomic_units(sysAU);

fprintf('\nAtomic-unit normalization:\n');
fprintf('  input length unit = angstrom\n');
fprintf('  output length unit = %s\n', sysAU.units.length);
fprintf('  output alpha unit  = %s\n', sysAU.units.alpha);
fprintf('  output charge unit = %s\n', sysAU.units.charge);

fprintf('\nFirst site position before conversion:\n');
disp(sys.site_pos(2,:));

fprintf('First site position after conversion:\n');
disp(sysAU.site_pos(2,:));

fprintf('\nTutorial 02 completed successfully.\n');

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end