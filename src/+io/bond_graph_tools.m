classdef bond_graph_tools
%BOND_GRAPH_TOOLS Shared low-level bond-graph utilities.
%
% Static methods:
%   r = io.bond_graph_tools.covalent_radius(sym)
%   tf = io.bond_graph_tools.ignore_pair(sym1, sym2)
%   cutoff = io.bond_graph_tools.bond_cutoff(sym1, sym2, scaleFactor)
%   tf = io.bond_graph_tools.is_bonded(sym1, sym2, d, scaleFactor)

    methods(Static)
        function r = covalent_radius(sym)
            switch upper(strtrim(sym))
                case 'H'
                    r = 0.31;
                case 'C'
                    r = 0.76;
                case 'N'
                    r = 0.71;
                case 'O'
                    r = 0.66;
                case 'F'
                    r = 0.57;
                case 'S'
                    r = 1.05;
                case 'CL'
                    r = 1.02;
                case 'BR'
                    r = 1.20;
                case 'I'
                    r = 1.39;
                case 'SI'
                    r = 1.11;
                otherwise
                    error('io:bond_graph_tools:UnknownElement', ...
                        'No covalent radius defined for element: %s', sym);
            end
        end

        function tf = ignore_pair(sym1, sym2)
        %IGNORE_PAIR Optional chemistry exclusions for bond finding.
            tf = strcmpi(strtrim(sym1), 'H') && strcmpi(strtrim(sym2), 'H');
        end

        function cutoff = bond_cutoff(sym1, sym2, scaleFactor)
            if nargin < 3 || isempty(scaleFactor)
                scaleFactor = 1.20;
            end

            if ~isscalar(scaleFactor) || ~isfinite(scaleFactor) || scaleFactor <= 0
                error('io:bond_graph_tools:BadScaleFactor', ...
                    'scaleFactor must be a positive finite scalar.');
            end

            r1 = io.bond_graph_tools.covalent_radius(sym1);
            r2 = io.bond_graph_tools.covalent_radius(sym2);
            cutoff = scaleFactor * (r1 + r2);
        end

        function tf = is_bonded(sym1, sym2, d, scaleFactor)
            if nargin < 4 || isempty(scaleFactor)
                scaleFactor = 1.20;
            end

            if ~isscalar(d) || ~isfinite(d) || d < 0
                error('io:bond_graph_tools:BadDistance', ...
                    'Distance d must be a nonnegative finite scalar.');
            end

            if io.bond_graph_tools.ignore_pair(sym1, sym2)
                tf = false;
                return;
            end

            cutoff = io.bond_graph_tools.bond_cutoff(sym1, sym2, scaleFactor);
            tf = (d < cutoff);
        end
    end
end