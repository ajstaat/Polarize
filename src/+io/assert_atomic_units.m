function assert_atomic_units(sys)
%ASSERT_ATOMIC_UNITS Ensure sys is normalized to internal atomic units.

    if ~isfield(sys, 'units') || isempty(sys.units)
        error('sys.units is missing; cannot verify unit convention.');
    end

    if ~isfield(sys.units, 'length') || ~strcmpi(sys.units.length, 'bohr')
        error('sys.site_pos must be in bohr.');
    end

    if ~isfield(sys.units, 'alpha') || ~strcmpi(sys.units.alpha, 'atomic_unit')
        error('sys.site_alpha must be in atomic units.');
    end

    if ~isfield(sys.units, 'charge') || ~strcmpi(sys.units.charge, 'elementary_charge')
        error('sys.site_charge must be in elementary-charge units.');
    end
end