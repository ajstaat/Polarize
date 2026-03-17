function validate_charge_pattern(chargePattern)
%VALIDATE_CHARGE_PATTERN Validate a site-label-based charge pattern.
%
% Required fields:
%   chargePattern.site_label   cell array of site labels
%   chargePattern.delta_q      numeric vector of same length
%
% Example:
%   chargePattern.site_label = {'C1','N1'};
%   chargePattern.delta_q = [0.12, -0.12];

if ~isstruct(chargePattern)
    error('chargePattern must be a struct.');
end

if ~isfield(chargePattern, 'site_label')
    error('chargePattern.site_label is required.');
end

if ~isfield(chargePattern, 'delta_q')
    error('chargePattern.delta_q is required.');
end

labels = chargePattern.site_label;
dq = chargePattern.delta_q;

if ~iscell(labels)
    error('chargePattern.site_label must be a cell array.');
end

if ~isnumeric(dq)
    error('chargePattern.delta_q must be numeric.');
end

if numel(labels) ~= numel(dq)
    error('chargePattern.site_label and chargePattern.delta_q must have the same length.');
end

if size(dq, 1) ~= 1 && size(dq, 2) ~= 1
    error('chargePattern.delta_q must be a vector.');
end

% Require unique labels within a pattern
if numel(unique(labels)) ~= numel(labels)
    error('chargePattern.site_label entries must be unique.');
end

end