function Tij = dipole_tensor_block(ri, rj, alpha_i, alpha_j, thole_a, opts)
%DIPOLE_TENSOR_BLOCK Dipole interaction tensor block from site j to site i.
%
% Inputs
%   ri, rj      1x3 Cartesian positions
%   alpha_i     scalar polarizability of site i
%   alpha_j     scalar polarizability of site j
%   thole_a     scalar Thole damping parameter
%   opts        struct with optional fields:
%               .softening   scalar, default 0
%               .use_thole   logical, default true
%
% Output
%   Tij         3x3 tensor such that E_i += Tij * mu_j

    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    softening = 0.0;
    if isfield(opts, 'softening') && ~isempty(opts.softening)
        softening = opts.softening;
    end

    use_thole = true;
    if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
        use_thole = opts.use_thole;
    end

    rvec = ri - rj;
    r2_bare = dot(rvec, rvec);

    if r2_bare == 0
        Tij = zeros(3,3);
        return;
    end

    r2 = r2_bare + softening^2;
    r  = sqrt(r2);
    r3 = r2 * r;
    r5 = r2 * r3;

    I3  = eye(3);
    rrT = rvec(:) * rvec(:).';

    if use_thole
        r_phys = sqrt(r2_bare);
        tf = thole.thole_f3f5_factors(r_phys, alpha_i, alpha_j, thole_a);
        Tij = (3 * tf.f5 / r5) * rrT - (tf.f3 / r3) * I3;
    else
        Tij = 3 * rrT / r5 - I3 / r3;
    end
end