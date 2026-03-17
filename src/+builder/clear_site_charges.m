function sys = clear_site_charges(sys)
%CLEAR_SITE_CHARGES Reset all site charges to zero.

if ~isfield(sys, 'n_sites') || isempty(sys.n_sites)
    error('sys.n_sites is missing or empty.');
end

sys.site_charge = zeros(sys.n_sites, 1);

end