function test_geom_spatial_index_nonperiodic()
%TEST_GEOM_SPATIAL_INDEX_NONPERIODIC Verify nonperiodic spatial-index queries.
%
% This compares brute-force and cell-list backends on a small nonperiodic
% Cartesian point cloud.

pos = [
    0.0 0.0 0.0
    1.0 0.0 0.0
    3.0 0.0 0.0
    0.0 2.0 0.0
    0.0 0.0 4.0
];

cutoff = 1.5;

opts = struct();
opts.isPeriodic = false;
opts.cutoff = cutoff;
opts.method = 'bruteforce';

spBF = geom.build_spatial_index(pos, opts);
pairsBF = geom.query_pairs_within_cutoff(spBF, cutoff, ...
    struct('return_r', true, 'return_dr', true));

opts.method = 'cell_list';

spCL = geom.build_spatial_index(pos, opts);
pairsCL = geom.query_pairs_within_cutoff(spCL, cutoff, ...
    struct('return_r', true, 'return_dr', true));

expectedPairs = [
    1 2
];

local_assert_pair_set(pairsBF.i, pairsBF.j, expectedPairs, ...
    'Brute-force nonperiodic pair set is wrong.');

local_assert_pair_set(pairsCL.i, pairsCL.j, expectedPairs, ...
    'Cell-list nonperiodic pair set is wrong.');

local_assert_same_pair_set(pairsBF, pairsCL, ...
    'Brute-force and cell-list nonperiodic pair sets differ.');

assert(numel(pairsBF.r) == 1 && abs(pairsBF.r(1) - 1.0) < 1e-12, ...
    'Expected brute-force nonperiodic pair distance is 1.0.');

assert(numel(pairsCL.r) == 1 && abs(pairsCL.r(1) - 1.0) < 1e-12, ...
    'Expected cell-list nonperiodic pair distance is 1.0.');

assert(isequal(size(pairsBF.dr), [1 3]), ...
    'Brute-force displacement array should be 1 x 3.');

assert(isequal(size(pairsCL.dr), [1 3]), ...
    'Cell-list displacement array should be 1 x 3.');

assert(norm(pairsBF.dr - [1.0 0.0 0.0], 'fro') < 1e-12, ...
    'Unexpected brute-force nonperiodic displacement.');

assert(norm(pairsCL.dr - [1.0 0.0 0.0], 'fro') < 1e-12, ...
    'Unexpected cell-list nonperiodic displacement.');

% Subset query: full indices [1,2,4] contain pairs (1,2) and (1,4)
% under cutoff 2.1.
subset = [1; 2; 4];
cutoff2 = 2.1;

opts.method = 'cell_list';
opts.cutoff = cutoff2;

sp = geom.build_spatial_index(pos, opts);

pairsSubsetLocal = geom.query_pairs_within_cutoff(sp, cutoff2, ...
    struct('subset_idx', subset, 'return_full_idx', false));

expectedLocal = [
    1 2
    1 3
];

local_assert_pair_set(pairsSubsetLocal.i, pairsSubsetLocal.j, expectedLocal, ...
    'Subset-local nonperiodic pair set is wrong.');

pairsSubsetFull = geom.query_pairs_within_cutoff(sp, cutoff2, ...
    struct('subset_idx', subset, 'return_full_idx', true));

expectedFull = [
    1 2
    1 4
];

local_assert_pair_set(pairsSubsetFull.i, pairsSubsetFull.j, expectedFull, ...
    'Subset-full nonperiodic pair set is wrong.');

end

function local_assert_same_pair_set(pairsA, pairsB, message)
A = sortrows([pairsA.i(:), pairsA.j(:)]);
B = sortrows([pairsB.i(:), pairsB.j(:)]);

assert(isequal(A, B), message);
end

function local_assert_pair_set(i, j, expected, message)
actual = sortrows([i(:), j(:)]);
expected = sortrows(expected);

assert(isequal(actual, expected), message);
end