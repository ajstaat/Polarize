function test_io_read_vasp_structure()
%TEST_IO_READ_VASP_STRUCTURE Verify minimal VASP POSCAR/CONTCAR parsing.
%
% This test creates temporary VASP-style files so the parser can be tested
% without relying on external crystal files.

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() local_cleanup(tmpDir));

%% Direct-coordinate POSCAR

directFile = fullfile(tmpDir, 'POSCAR_direct');

fid = fopen(directFile, 'w');
assert(fid > 0, 'Failed to open temporary direct POSCAR file.');

fprintf(fid, 'Tiny direct test\n');
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

S = io.read_vasp_structure(directFile);

assert(strcmp(S.comment, 'Tiny direct test'), ...
    'Comment line was not parsed correctly.');

assert(abs(S.scale - 1.0) < 1e-12, ...
    'Scale factor was not parsed correctly.');

assert(isequal(size(S.lattice), [3 3]), ...
    'Lattice should be 3 x 3.');

assert(S.natoms == 3, ...
    'Incorrect total atom count.');

assert(strcmp(S.coord_type, 'direct'), ...
    'Coordinate type should be direct.');

assert(numel(S.species) == 3, ...
    'Species label count should match atom count.');

assert(strcmp(S.species{1}, 'C'), ...
    'First species label should be C.');

assert(strcmp(S.species{2}, 'H') && strcmp(S.species{3}, 'H'), ...
    'Second and third species labels should be H.');

assert(norm(S.frac - [
    0.0 0.0 0.0
    0.5 0.0 0.0
    0.0 0.5 0.0
], 'fro') < 1e-12, ...
    'Direct fractional coordinates were parsed incorrectly.');

expectedCart = S.frac * S.lattice;

assert(norm(S.cart - expectedCart, 'fro') < 1e-12, ...
    'Direct-coordinate Cartesian conversion is incorrect.');

%% Cartesian-coordinate POSCAR

cartFile = fullfile(tmpDir, 'POSCAR_cart');

fid = fopen(cartFile, 'w');
assert(fid > 0, 'Failed to open temporary Cartesian POSCAR file.');

fprintf(fid, 'Tiny cartesian test\n');
fprintf(fid, '1.0\n');
fprintf(fid, '10.0 0.0 0.0\n');
fprintf(fid, '0.0 10.0 0.0\n');
fprintf(fid, '0.0 0.0 10.0\n');
fprintf(fid, 'O\n');
fprintf(fid, '2\n');
fprintf(fid, 'Cartesian\n');
fprintf(fid, '1.0 2.0 3.0\n');
fprintf(fid, '4.0 5.0 6.0\n');
fclose(fid);

S = io.read_vasp_structure(cartFile);

assert(S.natoms == 2, ...
    'Incorrect atom count for Cartesian POSCAR.');

assert(strcmp(S.coord_type, 'cartesian'), ...
    'Coordinate type should be cartesian.');

expectedCart = [
    1.0 2.0 3.0
    4.0 5.0 6.0
];

assert(norm(S.cart - expectedCart, 'fro') < 1e-12, ...
    'Cartesian coordinates were parsed incorrectly.');

expectedFrac = expectedCart / S.lattice;

assert(norm(S.frac - expectedFrac, 'fro') < 1e-12, ...
    'Cartesian-coordinate fractional conversion is incorrect.');

end

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end