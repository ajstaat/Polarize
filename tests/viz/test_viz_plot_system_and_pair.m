function test_viz_plot_system_and_pair()
%TEST_VIZ_PLOT_SYSTEM_AND_PAIR Smoke tests for system and molecule pair plotting.

sys = local_make_viz_sys();

fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>

ax = axes(fig);

out = viz.plot_system(sys, ...
    'Axes', ax, ...
    'ColorBy', 'molecule', ...
    'ShowBonds', true, ...
    'DrawBox', true, ...
    'Title', 'viz smoke test');

assert(isfield(out, 'ax') && isgraphics(out.ax), ...
    'plot_system should return a valid axes handle.');

assert(isfield(out, 'site_handle') && isgraphics(out.site_handle), ...
    'plot_system should return a site scatter handle.');

fig2 = figure('Visible', 'off');
cleanup2 = onCleanup(@() close(fig2)); %#ok<NASGU>

ax2 = axes(fig2);

pairOut = viz.plot_molecule_pair(sys, [1 2], ...
    'Axes', ax2, ...
    'RequireComplete', true);

assert(isfield(pairOut, 'ax') && isgraphics(pairOut.ax), ...
    'plot_molecule_pair should return a valid axes handle.');

assert(isequal(pairOut.mol_ids, [1 2]), ...
    'plot_molecule_pair should preserve molecule IDs.');

assert(numel(pairOut.site_indices) == 2, ...
    'plot_molecule_pair should return two site-index sets.');

end

function sys = local_make_viz_sys()

BOHR_PER_ANG = 1.8897259886;

% Two tiny ethane-ish C-H fragments far apart. Coordinates are in bohr.
X1A = [
    0.00 0.00 0.00
    1.09 0.00 0.00
] * BOHR_PER_ANG;

X2A = [
    4.00 0.00 0.00
    5.09 0.00 0.00
] * BOHR_PER_ANG;

site_pos = [X1A; X2A];

sys = struct();
sys.site_pos = site_pos;
sys.site_type = {'C'; 'H'; 'C'; 'H'};
sys.site_mol_id = [1; 1; 2; 2];
sys.unique_mol_id = sys.site_mol_id;

sys.units.length = 'bohr';
sys.super_lattice = 12 * BOHR_PER_ANG * eye(3);
sys.lattice = sys.super_lattice;

T = struct();
T.molecule_id = [1; 2];
T.unique_mol_id = [1; 2];
T.site_indices = {(1:2).'; (3:4).'};
T.n_sites = [2; 2];
T.com = [
    mean(site_pos(1:2,:), 1)
    mean(site_pos(3:4,:), 1)
];
T.is_complete_in_display = [true; true];

sys.molecule_table = T;
sys.site_charge = zeros(4,1);

end