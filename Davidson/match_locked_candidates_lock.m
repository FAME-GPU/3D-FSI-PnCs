function [is_locked_replica, lock_id, match_overlap] = match_locked_candidates_lock(cands, locked, opts)

    nc = numel(cands);
    nl = numel(locked.theta);
    is_locked_replica = false(nc,1);
    lock_id = zeros(nc,1);
    match_overlap = zeros(nc,1);

    if nc == 0 || nl == 0
        return;
    end

    thetaTol = get_lock_field_lock(opts,'thetaTol',1e-6);
    vectorTol = get_lock_field_lock(opts,'vectorTol',1e-3);
    min_overlap = 1-vectorTol;

    edges = zeros(0,3); 
    for j = 1:nc
        th = cands(j).theta;
        for l = 1:nl
            thL = locked.theta(l);
            scale = max([abs(th),abs(thL),eps]);
            theta_rel = abs(th-thL)/scale;
            if theta_rel > thetaTol
                continue;
            end

            ov = abs(locked.MsU_s(:,l)'*cands(j).u_s + ...
                     locked.MfU_f(:,l)'*cands(j).u_f);
            if isfinite(ov) && ov >= min_overlap
                edges(end+1,:) = [ov,j,l]; 
            end
        end
    end

    if isempty(edges)
        return;
    end

    [~,ord] = sort(edges(:,1),'descend');
    edges = edges(ord,:);
    used_c = false(nc,1);
    used_l = false(nl,1);

    for k = 1:size(edges,1)
        ov = edges(k,1);
        j = edges(k,2);
        l = edges(k,3);
        if ~used_c(j) && ~used_l(l)
            used_c(j) = true;
            used_l(l) = true;
            is_locked_replica(j) = true;
            lock_id(j) = l;
            match_overlap(j) = ov;
        end
    end
end

function value = get_lock_field_lock(opts,name,default_value)
if isfield(opts,'locking') && isstruct(opts.locking) && ...
        isfield(opts.locking,name) && ~isempty(opts.locking.(name))
    value = opts.locking.(name);
else
    value = default_value;
end
end
