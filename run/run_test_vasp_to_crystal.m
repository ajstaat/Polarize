clear; clc;

% ------------------------------------------------------------
% Build neutral system
% ------------------------------------------------------------
rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

molecules = io.unwrap_all_contcar_molecules(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', true);

viz.plot_all_unwrapped_contcar_molecules(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', true);

% clear; clc;
% 
% % ============================================================
% % Full test script:
% %   1) import CONTCAR
% %   2) build neutral supercell
% %   3) choose central cation/anion pair
% %   4) distribute -1/+1 uniformly over the two active molecules
% %   5) depolarize the active molecules
% %   6) run direct polarization calculation
% %   7) diagnose shortest charged->polarizable contact
% %   8) inspect suspicious molecule pair
% %   9) plot system / pair
% % ============================================================
% 
% %% -----------------------------
% % File / import
% % ------------------------------
% rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
% filename = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');
% 
% crystal = io.import_contcar_as_crystal(filename, ...
%     'BondScale', 1.20, ...
%     'WrapFractional', true, ...
%     'SortMolecules', false);
% 
% %% -----------------------------
% % Model
% % AMOEBA-like aromatic-ish test values
% % ------------------------------
% model = struct();
% model.polarizable_types = {'H','C','N','O'};
% model.alpha_by_type = struct( ...
%     'H', 0.696, ...
%     'C', 1.750, ...
%     'N', 1.073, ...
%     'O', 0.837);
% model.thole_a = 0.39;
% 
% % Neutral first pass
% model.charge_pattern = [];
% if isfield(model, 'charge_patterns')
%     model = rmfield(model, 'charge_patterns');
% end
% 
% %% -----------------------------
% % Options / params
% % ------------------------------
% opts = struct();
% opts.supercell_size = [2 2 1];
% opts.verbose = true;
% opts.removeMolIDs = [];
% opts.activeMolIDs = [];
% opts.depolarizeActiveMolecules = true;   % important
% 
% params = util.default_params();
% params.scf.solver = 'direct';            % use trusted reference solver
% params.scf.verbose = true;
% params.scf.printEvery = 10;
% params.scf.mixing = 0.05;
% params.field.softening = 0.5;
% params.scf.softening = 0.5;
% 
% %% -----------------------------
% % Neutral build
% % ------------------------------
% fprintf('\n--- Neutral build ---\n');
% result0 = calc.run_polarization_calc(crystal, model, opts, params);
% sys0 = result0.sys;
% 
% %% -----------------------------
% % Choose two molecules near center
% % ------------------------------
% [molCat, molAn, pairInfo] = builder.find_central_molecule_pair(sys0, ...
%     'DifferentBaseMol', true);
% 
% fprintf('Selected central pair: molCat=%d, molAn=%d\n', molCat, molAn);
% disp(pairInfo.rowA)
% disp(pairInfo.rowB)
% 
% result0 = calc.run_polarization_calc(crystal, model, opts, params);
% viz.plot_molecule_pair(result0.sys, 12, 16, ...
%     'HighlightChargedSites', false, ...
%     'ShowCOM', true, ...
%     'ShowMolIDs', true);
% 
% %% -----------------------------
% % Build uniform +/-1 charge patterns
% % ------------------------------
% idxCat = (sys0.site_mol_id == molCat);
% idxAn  = (sys0.site_mol_id == molAn);
% 
% labelsCat = sys0.site_label(idxCat);
% labelsAn  = sys0.site_label(idxAn);
% 
% nCat = numel(labelsCat);
% nAn  = numel(labelsAn);
% 
% fprintf('molCat has %d atoms, molAn has %d atoms\n', nCat, nAn);
% 
% % Cation/anion assignment:
% %   molCat gets -1 total
% %   molAn  gets +1 total
% model.charge_patterns = {
%     struct('mol_id', molCat, ...
%            'site_label', {labelsCat(:).'}, ...
%            'delta_q',    (-1/nCat) * ones(1, nCat)), ...
%     struct('mol_id', molAn, ...
%            'site_label', {labelsAn(:).'}, ...
%            'delta_q',    (+1/nAn) * ones(1, nAn))
% };
% 
% opts.activeMolIDs = [molCat, molAn];
% 
% %% -----------------------------
% % Charged run
% % ------------------------------
% fprintf('\n--- Charged run ---\n');
% result = calc.run_polarization_calc(crystal, model, opts, params);
% 
% %% -----------------------------
% % Energy summary
% % ------------------------------
% fprintf('\n--- Energy summary ---\n');
% disp(result.energy)
% 
% %% -----------------------------
% % Check active molecules really got depolarized
% % ------------------------------
% fprintf('\n--- Active-molecule depolarization check ---\n');
% idxActive = ismember(result.sys.site_mol_id, opts.activeMolIDs);
% 
% Tactive = table( ...
%     result.sys.site_mol_id(idxActive), ...
%     result.sys.site_label(idxActive), ...
%     result.sys.site_type(idxActive), ...
%     result.sys.site_charge(idxActive), ...
%     result.sys.site_alpha(idxActive), ...
%     result.sys.site_is_polarizable(idxActive), ...
%     'VariableNames', {'MolID','Label','Type','Charge','Alpha','IsPolarizable'});
% 
% disp(Tactive(1:min(20,height(Tactive)), :))
% 
% %% -----------------------------
% % Diagnose nearest charged -> polarizable contact
% % ------------------------------
% fprintf('\n--- Charged-to-polarizable contact diagnostic ---\n');
% 
% idxCharged = abs(result.sys.site_charge) > 1e-14;
% idxPol = logical(result.sys.site_is_polarizable(:));
% 
% xyzQ = result.sys.site_pos(idxCharged,:);
% xyzPol = result.sys.site_pos(idxPol,:);
% 
% allChargedIdx = find(idxCharged);
% allPolIdx = find(idxPol);
% 
% Dqp = pdist2(xyzQ, xyzPol);
% Dqp(Dqp < 1e-12) = inf;
% 
% [minQP, linIdx] = min(Dqp(:));
% [iq, ip] = ind2sub(size(Dqp), linIdx);
% 
% gq = allChargedIdx(iq);
% gp = allPolIdx(ip);
% 
% uidQ = result.sys.site_mol_id(gq);
% uidP = result.sys.site_mol_id(gp);
% 
% fprintf('Min charged-site to polarizable-site distance = %.4f A\n', minQP);
% fprintf('Charged site: %s (%s), mol %d\n', ...
%     result.sys.site_label{gq}, result.sys.site_type{gq}, uidQ);
% fprintf('Polarizable site: %s (%s), mol %d\n', ...
%     result.sys.site_label{gp}, result.sys.site_type{gp}, uidP);
% 
% disp(result.sys.site_pos(gq,:))
% disp(result.sys.site_pos(gp,:))
% 
% %% -----------------------------
% % Full molecule-pair minimum distance for suspicious pair
% % ------------------------------
% fprintf('\n--- Molecule-pair minimum-distance diagnostic ---\n');
% 
% idxQmol = (result.sys.site_mol_id == uidQ);
% idxPmol = (result.sys.site_mol_id == uidP);
% 
% xyzQmol = result.sys.site_pos(idxQmol,:);
% xyzPmol = result.sys.site_pos(idxPmol,:);
% 
% Dmol = pdist2(xyzQmol, xyzPmol);
% [minMol, linMol] = min(Dmol(:));
% [iqm, ipm] = ind2sub(size(Dmol), linMol);
% 
% globalQ = find(idxQmol);
% globalP = find(idxPmol);
% gqm = globalQ(iqm);
% gpm = globalP(ipm);
% 
% fprintf('Min atom-atom distance between mol %d and mol %d = %.4f A\n', ...
%     uidQ, uidP, minMol);
% fprintf('Site on mol %d: %s (%s)\n', uidQ, ...
%     result.sys.site_label{gqm}, result.sys.site_type{gqm});
% fprintf('Site on mol %d: %s (%s)\n', uidP, ...
%     result.sys.site_label{gpm}, result.sys.site_type{gpm});
% 
% disp(result.sys.site_pos(gqm,:))
% disp(result.sys.site_pos(gpm,:))
% 
% %% -----------------------------
% % COM table rows for suspicious pair
% % ------------------------------
% fprintf('\n--- Molecule table rows for suspicious pair ---\n');
% Tmol = builder.list_molecules(result.sys);
% disp(Tmol(Tmol.unique_mol_id == uidQ, :))
% disp(Tmol(Tmol.unique_mol_id == uidP, :))
% 
% %% -----------------------------
% % Optional visual checks
% % ------------------------------
% fprintf('\n--- Plotting ---\n');
% 
% viz.plot_molecule_pair(result.sys, molCat, molAn, ...
%     'ShowCOM', false, ...
%     'ShowMolIDs', false, ...
%     'HighlightChargedSites', false, ...
%     'ShowNeighbors', false, ...
%     'DrawCellBox', false, ...
%     'ShowLabels', false);
% 
% %% -----------------------------
% % Optional: save quick summary variables
% % ------------------------------
% summary = struct();
% summary.molCat = molCat;
% summary.molAn = molAn;
% summary.uidClosestChargedMol = uidQ;
% summary.uidClosestPolarizableMol = uidP;
% summary.minChargedToPolarizable = minQP;
% summary.minMolPairDistance = minMol;
% summary.energy = result.energy;
% 
% fprintf('\n--- Summary ---\n');
% disp(summary)