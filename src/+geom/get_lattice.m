function lat = get_lattice(sysOrHOrLat)
%GET_LATTICE Return canonical project lattice struct.
%
% Project convention:
%
%   H has direct lattice vectors as rows:
%
%       r_cart = f_frac * H
%
%   G has reciprocal lattice vectors as columns:
%
%       H * G = 2*pi*I
%
%   For integer reciprocal index m_col:
%
%       k_col = G * m_col
%       k_row = k_col.'
%
% Accepted inputs:
%
%   1. Existing lattice struct with fields:
%        .H
%        .G
%
%   2. System struct with either:
%        .super_lattice
%        .lattice
%
%   3. Raw 3x3 H matrix, direct lattice vectors as rows.
%
% Output:
%
%   lat.H
%   lat.G
%   lat.volume
%   lat.identity_error
%   lat.convention

    if nargin < 1 || isempty(sysOrHOrLat)
        error('geom:get_lattice:MissingInput', ...
            'Input must be a lattice struct, system struct, or 3x3 H matrix.');
    end

    % ---------------------------------------------------------------------
    % Case 1: already a canonical lattice struct.
    %
    % Important: allow this pass-through so downstream routines can safely do:
    %
    %   lat = geom.get_lattice(latOrSysOrH);
    %
    % even when latOrSysOrH is already lat.
    % ---------------------------------------------------------------------

    if isstruct(sysOrHOrLat) && ...
       isfield(sysOrHOrLat, 'H') && ~isempty(sysOrHOrLat.H) && ...
       isfield(sysOrHOrLat, 'G') && ~isempty(sysOrHOrLat.G)

        lat = sysOrHOrLat;

        validateattributes(lat.H, {'double'}, {'size',[3 3], 'real', 'finite'}, ...
            mfilename, 'lat.H');
        validateattributes(lat.G, {'double'}, {'size',[3 3], 'real', 'finite'}, ...
            mfilename, 'lat.G');

        if ~isfield(lat, 'volume') || isempty(lat.volume)
            lat.volume = abs(det(lat.H));
        end

        if ~isfield(lat, 'identity_error') || isempty(lat.identity_error)
            lat.identity_error = max(abs(lat.H * lat.G - 2*pi*eye(3)), [], 'all');
        end

        if ~isfield(lat, 'convention') || isempty(lat.convention)
            lat.convention = 'H_rows_G_columns_HG_2piI';
        end

        if lat.volume <= 1e-14
            error('geom:get_lattice:SingularCell', ...
                'Direct lattice has near-zero volume.');
        end

        if lat.identity_error > 1e-10
            error('geom:get_lattice:BadReciprocalIdentity', ...
                'Input lattice struct fails H*G = 2*pi*I; error %.3e.', ...
                lat.identity_error);
        end

        return;
    end

    % ---------------------------------------------------------------------
    % Case 2: system struct with .super_lattice or .lattice.
    % ---------------------------------------------------------------------

    if isstruct(sysOrHOrLat)
        sys = sysOrHOrLat;

        if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
            H = sys.super_lattice;
        elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
            H = sys.lattice;
        else
            error('geom:get_lattice:MissingLattice', ...
                'Missing sys.super_lattice or sys.lattice.');
        end

    % ---------------------------------------------------------------------
    % Case 3: raw 3x3 H matrix.
    % ---------------------------------------------------------------------

    else
        H = sysOrHOrLat;
    end

    validateattributes(H, {'double'}, {'size',[3 3], 'real', 'finite'}, ...
        mfilename, 'H');

    lat = struct();
    lat.H = H;
    lat.volume = abs(det(H));

    if lat.volume <= 1e-14
        error('geom:get_lattice:SingularCell', ...
            'Direct lattice has near-zero volume.');
    end

    % Project reciprocal convention:
    %
    %   H * G = 2*pi*I
    %
    % Therefore:
    %
    %   G = 2*pi * inv(H)
    %
    % Columns of G are reciprocal lattice vectors.
    lat.G = 2*pi * inv(H);

    lat.identity_error = max(abs(lat.H * lat.G - 2*pi*eye(3)), [], 'all');
    lat.convention = 'H_rows_G_columns_HG_2piI';

    if lat.identity_error > 1e-10
        error('geom:get_lattice:BadReciprocalIdentity', ...
            'H*G reciprocal identity failed: %.3e.', lat.identity_error);
    end
end