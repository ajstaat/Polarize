function Tfull = expand_active_operator_to_full(problem, Tpol)
    Tfull = zeros(3*problem.nSites, 3*problem.nSites);
    idx = problem.activeVecIdx;
    Tfull(idx, idx) = Tpol;
end