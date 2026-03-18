function pretty_print_result(result)
%PRETTY_PRINT_RESULT Print a concise summary of a calculation result.

fprintf('\nStatus: %s\n', result.status);
fprintf('%s\n', result.message);

if isfield(result, 'scf') && ~isempty(result.scf)
    if isfield(result.scf, 'converged')
        fprintf('SCF converged: %d\n', result.scf.converged);
    end
    if isfield(result.scf, 'nIter')
        fprintf('SCF iterations: %d\n', result.scf.nIter);
    end
end

if isfield(result, 'energy') && ~isempty(result.energy)
    e = result.energy;

    if isfield(e, 'polarization_self')
        fprintf('U_pol  = %+ .12e\n', e.polarization_self);
    end
    if isfield(e, 'external_charge_dipole')
        fprintf('U_ext  = %+ .12e\n', e.external_charge_dipole);
    end
    if isfield(e, 'dipole_dipole')
        fprintf('U_dd   = %+ .12e\n', e.dipole_dipole);
    end
    if isfield(e, 'total')
        fprintf('U_tot  = %+ .12e\n', e.total);
    end
end

if isfield(result, 'Tmeta') && ~isempty(result.Tmeta) && isstruct(result.Tmeta)
    if isfield(result.Tmeta, 'num_kvec')
        fprintf('num_kvec = %d\n', result.Tmeta.num_kvec);
    end
end

if isfield(result, 'dipole_dipole_decomposition') && ~isempty(result.dipole_dipole_decomposition)
    dd = result.dipole_dipole_decomposition;
    if isfield(dd, 'active_active')
        fprintf('U_dd(active,active)   = %+ .12e\n', dd.active_active);
        fprintf('U_dd(active,env)      = %+ .12e\n', dd.active_environment);
        fprintf('U_dd(env,env)         = %+ .12e\n', dd.environment_environment);
    end
end

end