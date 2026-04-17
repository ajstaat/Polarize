function S = read_vasp_structure(filename)
%READ_VASP_STRUCTURE Parse a VASP POSCAR/CONTCAR-style structure file.
%
% Returns fields:
%   comment
%   scale
%   lattice
%   species
%   natoms
%   frac
%   cart
%   coord_type   ('direct' or 'cartesian')

    fid = fopen(filename, 'r');
    if fid < 0
        error('io:read_vasp_structure:FileOpenFailed', ...
            'Could not open file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            lines{end+1,1} = strtrim(line); %#ok<AGROW>
        end
    end

    if numel(lines) < 8
        error('io:read_vasp_structure:FileTooShort', ...
            'File too short to be a valid VASP structure: %s', filename);
    end

    S = struct();
    S.comment = lines{1};

    S.scale = str2double(lines{2});
    if ~isfinite(S.scale) || S.scale == 0
        error('io:read_vasp_structure:BadScale', ...
            'Invalid scale factor in file: %s', filename);
    end

    lattice = zeros(3,3);
    for i = 1:3
        vals = sscanf(lines{2+i}, '%f %f %f');
        if numel(vals) < 3
            error('io:read_vasp_structure:BadLattice', ...
                'Failed to parse lattice vector %d in %s', i, filename);
        end
        lattice(i,:) = vals(1:3).';
    end
    S.lattice = S.scale * lattice;

    speciesLine = strsplit(strtrim(lines{6}));
    countsLine  = sscanf(lines{7}, '%d').';

    if isempty(speciesLine) || isempty(countsLine)
        error('io:read_vasp_structure:BadSpeciesCounts', ...
            'Could not parse species/counts in file: %s', filename);
    end
    if numel(speciesLine) ~= numel(countsLine)
        error('io:read_vasp_structure:SpeciesCountMismatch', ...
            'Species/count mismatch in file: %s', filename);
    end
    if any(countsLine < 0)
        error('io:read_vasp_structure:NegativeCounts', ...
            'Negative atom count found in file: %s', filename);
    end

    natoms = sum(countsLine);
    S.natoms = natoms;

    coordLineIdx = 8;
    if startsWith(lower(lines{8}), 's')
        coordLineIdx = 9;
    end

    if coordLineIdx > numel(lines)
        error('io:read_vasp_structure:MissingCoordinateType', ...
            'Missing coordinate type line in file: %s', filename);
    end

    coordType = lower(strtrim(lines{coordLineIdx}));
    if startsWith(coordType, 'd')
        isDirect = true;
        S.coord_type = 'direct';
    elseif startsWith(coordType, 'c')
        isDirect = false;
        S.coord_type = 'cartesian';
    else
        error('io:read_vasp_structure:UnknownCoordinateType', ...
            'Could not determine coordinate type in file: %s', filename);
    end

    startIdx = coordLineIdx + 1;
    endIdx   = startIdx + natoms - 1;
    if endIdx > numel(lines)
        error('io:read_vasp_structure:MissingCoordinates', ...
            'Not enough coordinate lines in file: %s', filename);
    end

    coords = zeros(natoms,3);
    for i = 1:natoms
        vals = sscanf(lines{startIdx+i-1}, '%f %f %f');
        if numel(vals) < 3
            error('io:read_vasp_structure:BadCoordinateLine', ...
                'Bad coordinate line %d in %s', i, filename);
        end
        coords(i,:) = vals(1:3).';
    end

    labels = cell(natoms,1);
    c = 1;
    for s = 1:numel(speciesLine)
        for j = 1:countsLine(s)
            labels{c} = speciesLine{s};
            c = c + 1;
        end
    end

    if isDirect
        frac = coords;
        cart = frac * S.lattice;
    else
        cart = coords;
        frac = cart / S.lattice;
    end

    S.species = labels;
    S.frac    = frac;
    S.cart    = cart;
end