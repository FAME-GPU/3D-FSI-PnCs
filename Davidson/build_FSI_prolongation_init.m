function [Vs0, Vf0, info] = build_FSI_prolongation_init( ...
    mesh_c, mesh_h, T_s_c, T_f_c, T_s_h, T_f_h, Phi_s_c, Phi_f_c, M_s_h, M_f_h, opts)

    if nargin < 11 || isempty(opts)
        opts = struct();
    end
    opts = set_default_prolong_opts(opts);

    if ~isfield(opts, 'computeSolid') || isempty(opts.computeSolid)
        opts.computeSolid = ~isempty(Phi_s_c);
    end
    if ~isfield(opts, 'computeFluid') || isempty(opts.computeFluid)
        opts.computeFluid = ~isempty(Phi_f_c);
    end

    Vs0 = [];
    Vf0 = [];

    info = struct();
    info.opts = opts;
    info.solid = struct('skipped', ~opts.computeSolid);
    info.fluid = struct('skipped', ~opts.computeFluid);
    info.size = struct();

    if ~opts.computeSolid && ~opts.computeFluid
        if opts.verbose
            fprintf('[Prolong] both solid and fluid prolongation are skipped. Return Vs0=[], Vf0=[].\n');
        end
        info.size.Vs0 = size(Vs0);
        info.size.Vf0 = size(Vf0);
        return;
    end

    if opts.computeSolid
        validate_part_inputs('solid', Phi_s_c, M_s_h, T_s_c, T_s_h, opts);

        if opts.verbose
            fprintf('[Prolong] building solid P1 prolongation...\n');
        end

        [coord_s_c, tet_s_c, coord_s_h] = extract_domain_geometry(mesh_c, mesh_h, 'solid');
        [P_s_node, info_s] = build_P1_prolongation(coord_s_c, tet_s_c, coord_s_h, opts);
        P_s_tight = kron(P_s_node, speye(3));  

        switch lower(opts.inputSpace)
            case 'reduced'
                Phi_s_c_tight = T_s_c * Phi_s_c;
            case 'tight'
                Phi_s_c_tight = Phi_s_c;
            otherwise
                error('opts.inputSpace must be ''reduced'' or ''tight''.');
        end

        Phi_s_h_tight = P_s_tight * Phi_s_c_tight;
        Vs0 = project_tight_to_reduced(T_s_h, Phi_s_h_tight);

        if opts.massOrth
            Vs0 = m_orthonormalize(Vs0, M_s_h, opts.orthTol, opts.maxReorth);
        end

        info.solid = info_s;
        info.solid.skipped = false;
        info.size.P_s_node = size(P_s_node);
        info.size.P_s_tight = size(P_s_tight);
        info.size.Vs0 = size(Vs0);
    else
        if opts.verbose
            fprintf('[Prolong] solid prolongation skipped. Vs0 = [].\n');
        end
        info.size.Vs0 = size(Vs0);
    end

    if opts.computeFluid
        validate_part_inputs('fluid', Phi_f_c, M_f_h, T_f_c, T_f_h, opts);

        if opts.verbose
            fprintf('[Prolong] building fluid P1 prolongation...\n');
        end

        [coord_f_c, tet_f_c, coord_f_h] = extract_domain_geometry(mesh_c, mesh_h, 'fluid');
        [P_f_tight, info_f] = build_P1_prolongation(coord_f_c, tet_f_c, coord_f_h, opts);

        switch lower(opts.inputSpace)
            case 'reduced'
                Phi_f_c_tight = T_f_c * Phi_f_c;
            case 'tight'
                Phi_f_c_tight = Phi_f_c;
            otherwise
                error('opts.inputSpace must be ''reduced'' or ''tight''.');
        end

        Phi_f_h_tight = P_f_tight * Phi_f_c_tight;
        Vf0 = project_tight_to_reduced(T_f_h, Phi_f_h_tight);

        if opts.massOrth
            Vf0 = m_orthonormalize(Vf0, M_f_h, opts.orthTol, opts.maxReorth);
        end

        info.fluid = info_f;
        info.fluid.skipped = false;
        info.size.P_f_tight = size(P_f_tight);
        info.size.Vf0 = size(Vf0);
    else
        if opts.verbose
            fprintf('[Prolong] fluid prolongation skipped. Vf0 = [].\n');
        end
        info.size.Vf0 = size(Vf0);
    end

    if opts.verbose
        fprintf('[Prolong] done. Vs0 dim = %d, Vf0 dim = %d.\n', size(Vs0,2), size(Vf0,2));
    end
end

function opts = set_default_prolong_opts(opts)
    opts = setfield_if_missing(opts, 'inputSpace', 'reduced');
    opts = setfield_if_missing(opts, 'massOrth', true);
    opts = setfield_if_missing(opts, 'orthTol', 1e-12);
    opts = setfield_if_missing(opts, 'maxReorth', 5);
    opts = setfield_if_missing(opts, 'locateTol', 1e-10);
    opts = setfield_if_missing(opts, 'maxFallbackBad', 2000);
    opts = setfield_if_missing(opts, 'fallbackTets', 80);
    opts = setfield_if_missing(opts, 'verbose', true);
end

function s = setfield_if_missing(s, name, val)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = val;
    end
end

function validate_part_inputs(domain, Phi_c, M_h, T_c, T_h, opts)

    if isempty(Phi_c)
        error('opts.compute%s=true, but Phi_%s_c is empty.', capitalize(domain), domain(1));
    end

    if opts.massOrth && isempty(M_h)
        error('opts.compute%s=true and opts.massOrth=true, but M_%s_h is empty.', ...
            capitalize(domain), domain(1));
    end

    switch lower(opts.inputSpace)
        case 'reduced'
            if isempty(T_c) || isempty(T_h)
                error(['opts.inputSpace=''reduced'' requires both coarse and fine Bloch maps ', ...
                       'for the %s part.'], domain);
            end
        case 'tight'
            if isempty(T_h)
                error('opts.inputSpace=''tight'' still requires the fine Bloch map T_%s_h.', domain(1));
            end
        otherwise
            error('opts.inputSpace must be ''reduced'' or ''tight''.');
    end
end

function out = capitalize(str)
    out = str;
    if ~isempty(out)
        out(1) = upper(out(1));
    end
end

function [coord_c, tet_c_local, coord_h] = extract_domain_geometry(mesh_c, mesh_h, domain)

    switch lower(domain)
        case 'solid'
            elem_c_idx = find(mesh_c.cellIndices == mesh_c.solidTag);
            nodes_c = mesh_c.nodeidx_s(:);
            nodes_h = mesh_h.nodeidx_s(:);
        case 'fluid'
            elem_c_idx = find(mesh_c.cellIndices == mesh_c.fluidTag);
            nodes_c = mesh_c.nodeidx_f(:);
            nodes_h = mesh_h.nodeidx_f(:);
        otherwise
            error('domain must be solid or fluid.');
    end

    tet_c_global = mesh_c.elements(:, elem_c_idx).';   

    maxNode = max([nodes_c(:); tet_c_global(:)]);
    g2l = zeros(maxNode, 1);
    g2l(nodes_c) = 1:numel(nodes_c);
    tet_c_local = g2l(tet_c_global);

    if any(tet_c_local(:) == 0)
        error('Some tetrahedron nodes are not contained in meshData.nodeidx_%s.', domain(1));
    end

    coord_c = mesh_c.nodes(:, nodes_c).';  
    coord_h = mesh_h.nodes(:, nodes_h).';  
end

function [P, info] = build_P1_prolongation(coord_c, tet_c, coord_h, opts)

    coord_c = double(coord_c);
    coord_h = double(coord_h);
    tet_c = double(tet_c);

    TR = triangulation(tet_c, coord_c);
    elem_id = pointLocation(TR, coord_h);

    bad = isnan(elem_id);
    n_bad_initial = nnz(bad);

    info = struct();
    info.nFine = size(coord_h, 1);
    info.nCoarse = size(coord_c, 1);
    info.nTetCoarse = size(tet_c, 1);
    info.nBadInitial = n_bad_initial;
    info.nBadAfterFallback = 0;

    bary = nan(size(coord_h,1), 4);
    good = ~bad;
    if any(good)
        bary(good, :) = cartesianToBarycentric(TR, elem_id(good), coord_h(good,:));
    end

    if n_bad_initial > 0
        if opts.verbose
            fprintf('[Prolong] pointLocation failed for %d/%d nodes; trying fallback.\n', ...
                n_bad_initial, size(coord_h,1));
        end
        if n_bad_initial > opts.maxFallbackBad
            error(['Too many fine nodes were not located in the coarse mesh. ', ...
                   'Use a nested refined mesh or increase opts.maxFallbackBad intentionally.']);
        end

        cent = (coord_c(tet_c(:,1),:) + coord_c(tet_c(:,2),:) + ...
                coord_c(tet_c(:,3),:) + coord_c(tet_c(:,4),:)) / 4;
        bad_idx = find(bad);
        for kk = 1:numel(bad_idx)
            ip = bad_idx(kk);
            [eid, b] = fallback_locate_one_point(coord_h(ip,:), coord_c, tet_c, cent, opts);
            if eid > 0
                elem_id(ip) = eid;
                bary(ip,:) = b;
                bad(ip) = false;
            end
        end
    end

    if any(bad)
        info.nBadAfterFallback = nnz(bad);
        bad_idx = find(bad);
        error('P1 prolongation failed: %d fine nodes not found. First bad local node = %d.', ...
            nnz(bad), bad_idx(1));
    end

    tol = opts.locateTol;
    bary(abs(bary) < tol) = 0;
    bary(abs(bary - 1) < tol) = 1;
    rowsum = sum(bary, 2);
    bary = bary ./ rowsum;

    n_h = size(coord_h, 1);
    n_c = size(coord_c, 1);

    I = repelem((1:n_h)', 4);
    J = reshape(tet_c(elem_id, :)', [], 1);
    S = reshape(bary', [], 1);

    P = sparse(I, J, S, n_h, n_c);
    info.nBadAfterFallback = nnz(bad);
end

function [eid, b] = fallback_locate_one_point(p, coord, tet, cent, opts)

    d2 = sum((cent - p).^2, 2);
    ktry = min(opts.fallbackTets, numel(d2));

    if exist('mink', 'builtin') || exist('mink', 'file')
        [~, cand] = mink(d2, ktry);
    else
        [~, ord] = sort(d2, 'ascend');
        cand = ord(1:ktry);
    end

    tol = opts.locateTol;
    eid = 0;
    b = [];

    for jj = 1:numel(cand)
        e = cand(jj);
        X = coord(tet(e,:), :);  
        A = [X(2,:) - X(1,:); X(3,:) - X(1,:); X(4,:) - X(1,:)].';
        rhs = (p - X(1,:)).';
        if rcond(A) < 1e-14
            continue;
        end
        xi = A \ rhs;
        bt = [1 - sum(xi), xi.'];
        if all(bt >= -tol) && all(bt <= 1 + tol)
            bt(abs(bt) < tol) = 0;
            bt(abs(bt - 1) < tol) = 1;
            b = bt / sum(bt);
            eid = e;
            return;
        end
    end
end

function Xred = project_tight_to_reduced(T, Xtight)

    d = full(sum(abs(T).^2, 1)).';
    if all(d > 0)
        Xred = T' * Xtight;
        Xred = bsxfun(@rdivide, Xred, d);
    else
        Xred = (T' * T) \ (T' * Xtight);
    end
end

function Q = m_orthonormalize(X, M, tol, maxReorth)

    if nargin < 3 || isempty(tol)
        tol = 1e-12;
    end
    if nargin < 4 || isempty(maxReorth)
        maxReorth = 5;
    end

    if isempty(X)
        Q = X;
        return;
    end

    Q = zeros(size(X,1), 0, 'like', X);

    for j = 1:size(X,2)
        v = X(:,j);

        for rep = 1:maxReorth
            if ~isempty(Q)
                h = Q' * (M * v);
                v = v - Q * h;
            end
        end

        Mv = M * v;
        nv = sqrt(abs(v' * Mv));
        if isfinite(nv) && nv > tol
            Q(:,end+1) = v / nv; 
        end
    end
end
