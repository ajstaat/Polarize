function draw_cell_box(ax, lattice, origin, vectors_as)
%DRAW_CELL_BOX Draw a parallelepiped from lattice vectors.
%
% draw_cell_box(ax, lattice)
% draw_cell_box(ax, lattice, origin)
% draw_cell_box(ax, lattice, origin, vectors_as)
%
% vectors_as:
%   'rows'    -> lattice(1,:), lattice(2,:), lattice(3,:)
%   'columns' -> lattice(:,1), lattice(:,2), lattice(:,3)
%
% default = 'rows'

    if nargin < 3 || isempty(origin)
        origin = [0 0 0];
    end
    if nargin < 4 || isempty(vectors_as)
        vectors_as = 'rows';
    end

    switch lower(vectors_as)
        case 'rows'
            a = lattice(1,:);
            b = lattice(2,:);
            c = lattice(3,:);
        case 'columns'
            a = lattice(:,1).';
            b = lattice(:,2).';
            c = lattice(:,3).';
        otherwise
            error('vectors_as must be ''rows'' or ''columns''.');
    end

    p000 = origin;
    p100 = origin + a;
    p010 = origin + b;
    p001 = origin + c;
    p110 = origin + a + b;
    p101 = origin + a + c;
    p011 = origin + b + c;
    p111 = origin + a + b + c;

    edges = [
        p000; p100; p110; p010; p000; ...
        p001; p101; p100; ...
        p101; p111; p110; ...
        p111; p011; p010; ...
        p011; p001
    ];

    plot3(ax, edges(:,1), edges(:,2), edges(:,3), 'k-', 'LineWidth', 1.2);
end