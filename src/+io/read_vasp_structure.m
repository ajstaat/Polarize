% src/+io/read_vasp_structure.m
function S = read_vasp_structure(filename)
%READ_VASP_STRUCTURE Parse a VASP POSCAR/CONTCAR-style structure file.

    fid = fopen(filename, 'r');
    if fid < 0
        error('Could not open file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = {};
    while ~feof(fid)
        lines{end+1,1} = strtrim(fgetl(fid)); %#ok<AGROW>
    end

    if numel(lines) < 8
        error('File too short to be a valid VASP structure: %s', filename);
    end

    S = struct();
    S.comment = lines{1};
    S.scale   = str2double(lines{2});

    lattice = zeros(3,3);
    for i = 1:3
        vals = sscanf(lines{2+i}, '%f %f %f');
        if numel(vals) < 3
            error('Failed to parse lattice vector %d in %s', i, filename);
        end
        lattice(i,:) = vals(1:3).';
    end
    S.lattice = S.scale * lattice;

    speciesLine = strsplit(strtrim(lines{6}));
    countsLine  = sscanf(lines{7}, '%f').';

    if isempty(speciesLine) || isempty(countsLine)
        error('Could not parse species/counts in file: %s', filename);
    end
    if numel(speciesLine) ~= numel(countsLine)
        error('Species/count mismatch in file: %s', filename);
    end

    natoms = sum(countsLine);

    coordLineIdx = 8;
    if startsWith(lower(lines{8}), 's')
        coordLineIdx = 9;
    end

    coordType = lower(strtrim(lines{coordLineIdx}));
    if startsWith(coordType, 'd')
        isDirect = true;
    elseif startsWith(coordType, 'c')
        isDirect = false;
    else
        error('Could not determine coordinate type in file: %s', filename);
    end

    startIdx = coordLineIdx + 1;
    endIdx   = startIdx + natoms - 1;
    if endIdx > numel(lines)
        error('Not enough coordinate lines in file: %s', filename);
    end

    coords = zeros(natoms,3);
    for i = 1:natoms
        vals = sscanf(lines{startIdx+i-1}, '%f %f %f');
        if numel(vals) < 3
            error('Bad coordinate line %d in %s', i, filename);
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