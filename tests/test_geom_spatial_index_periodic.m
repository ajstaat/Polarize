function test_geom_spatial_index_periodic()
%TEST_GEOM_SPATIAL_INDEX_PERIODIC Verify periodic spatial-index queries.
%
% This test checks the row-vector lattice convention:
%
%   cart = frac * H
%
% It compares brute-force and cell-list backends for a boundary-crossing
% pair under minimum-image periodic geometry.

H = [
    10.0 0.0 0.0
     0.0 10.0 0.0
     0.0 0.0 10.0
];

pos = [
    9.5 5.0 5.0
    0.5 5.0 5.0
    5.0 5.0 5.0
    5.0 5.0 6.0
];

% Under PBC:
%   pair (1,2) has minimum-image distance 1.0 across the x boundary.
%   pair (3,4) has distance 1.0 inside the cell.
% Other pairs are farther than cutoff.
cutoff = 1.5;

opts = struct();
opts.isPeriodic = true;
opts.cell = H;
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
    3 4
];

local_assert_pair_set(pairsBF.i, pairsBF.j, expectedPairs, ...
    'Brute-force periodic pair set is wrong.');

local_assert_pair_set(pairsCL.i, pairsCL.j, expectedPairs, ...
    'Cell-list periodic pair set is wrong.');

local_assert_same_pair_set(pairsBF, pairsCL, ...
    'Brute-force and cell-list periodic pair sets differ.');

% Check distances independent of ordering.
local_assert_distance_for_pair(pairsBF, 1, 2, 1.0, ...
    'Brute-force periodic boundary distance should be 1.0.');

local_assert_distance_for_pair(pairsCL, 1, 2, 1.0, ...
    'Cell-list periodic boundary distance should be 1.0.');

local_assert_distance_for_pair(pairsBF, 3, 4, 1.0, ...
    'Brute-force internal periodic distance should be 1.0.');

local_assert_distance_for_pair(pairsCL, 3, 4, 1.0, ...
    'Cell-list internal periodic distance should be 1.0.');

% Check the boundary-crossing displacement is minimum-image. The pair is
% canonicalized as (1,2), so displacement should be from point 1 to point 2.
% From x=9.5 to x=0.5 under PBC, the minimum image is +1.0 along x.
local_assert_dr_for_pair(pairsBF, 1, 2, [1.0 0.0 0.0], ...
    'Unexpected brute-force periodic boundary displacement.');

local_assert_dr_for_pair(pairsCL, 1, 2, [1.0 0.0 0.0], ...
    'Unexpected cell-list periodic boundary displacement.');

% Also check a sheared row-lattice case. This catches row/column convention
% mistakes that cubic lattices may hide.
Hshear = [
    10.0  0.0  0.0
     2.0 10.0  0.0
     0.0  0.0 10.0
];

frac = [
    0.95 0.50 0.50
    0.05 0.50 0.50
    0.50 0.50 0.50
];

posShear = frac * Hshear;

opts = struct();
opts.isPeriodic = true;
opts.cell = Hshear;
opts.cutoff = cutoff;
opts.method = 'bruteforce';

spBF = geom.build_spatial_index(posShear, opts);
pairsBF = geom.query_pairs_within_cutoff(spBF, cutoff, ...
    struct('return_r', true, 'return_dr', true));

opts.method = 'cell_list';

spCL = geom.build_spatial_index(posShear, opts);
pairsCL = geom.query_pairs_within_cutoff(spCL, cutoff, ...
    struct('return_r', true, 'return_dr', true));

expectedPairs = [
    1 2
];

local_assert_pair_set(pairsBF.i, pairsBF.j, expectedPairs, ...
    'Brute-force sheared periodic pair set is wrong.');

local_assert_pair_set(pairsCL.i, pairsCL.j, expectedPairs, ...
    'Cell-list sheared periodic pair set is wrong.');

local_assert_same_pair_set(pairsBF, pairsCL, ...
    'Brute-force and cell-list sheared periodic pair sets differ.');

local_assert_distance_for_pair(pairsBF, 1, 2, 1.0, ...
    'Brute-force sheared periodic distance should be 1.0.');

local_assert_distance_for_pair(pairsCL, 1, 2, 1.0, ...
    'Cell-list sheared periodic distance should be 1.0.');

local_assert_dr_for_pair(pairsBF, 1, 2, [1.0 0.0 0.0], ...
    'Unexpected brute-force sheared periodic displacement.');

local_assert_dr_for_pair(pairsCL, 1, 2, [1.0 0.0 0.0], ...
    'Unexpected cell-list sheared periodic displacement.');

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

function local_assert_distance_for_pair(pairs, i, j, expectedDistance, message)
idx = local_find_pair(pairs, i, j);

assert(~isempty(idx), ...
    'Requested pair was not found.');

assert(abs(pairs.r(idx) - expectedDistance) < 1e-12, message);
end

function local_assert_dr_for_pair(pairs, i, j, expectedDr, message)
idx = local_find_pair(pairs, i, j);

assert(~isempty(idx), ...
    'Requested pair was not found.');

assert(norm(pairs.dr(idx,:) - expectedDr, 'fro') < 1e-12, message);
end

function idx = local_find_pair(pairs, i, j)
idx = find(pairs.i == i & pairs.j == j, 1, 'first');
end