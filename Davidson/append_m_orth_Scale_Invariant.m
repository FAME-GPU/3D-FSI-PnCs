function [V, added, diaginfo] = append_m_orth_Scale_Invariant(V, Z, M, opts)

    empty = struct('input_col',[],'finite',false,'base_norm',NaN, ...
                   'final_norm',NaN,'norm_ratio',NaN,'relative_tol',opts.orthTol, ...
                   'reorth_steps',0,'dgks_accept',false,'appended',false, ...
                   'reason','');
    diaginfo = repmat(empty,size(Z,2),1);

    if isempty(Z)
        added = 0;
        return;
    end

    added = 0;
    for j = 1:size(Z,2)
        diaginfo(j).input_col = j;
        z = Z(:,j);
        if any(~isfinite(z))
            diaginfo(j).reason = 'nonfinite_input';
            continue;
        end
        diaginfo(j).finite = true;

        Mz = M*z;
        base = sqrt(max(real(z'*Mz),0));
        diaginfo(j).base_norm = base;

        if ~(isfinite(base) && base > 0)
            diaginfo(j).final_norm = base;
            diaginfo(j).norm_ratio = NaN;
            diaginfo(j).reason = 'zero_or_invalid_base';
            continue;
        end

        norm_old = base;
        norm_new = base;
        accepted = false;

        for k = 1:opts.maxReorth
            if ~isempty(V)
                h = V'*Mz;
                z = z - V*h;
                Mz = M*z;
            end
            norm_new = sqrt(max(real(z'*Mz),0));
            diaginfo(j).reorth_steps = k;

            if norm_new > opts.orthTol*base && ...
               norm_new > (1/sqrt(2))*norm_old
                accepted = true;
                break;
            end
            norm_old = max(norm_new,realmin(class(real(norm_new))));
        end

        diaginfo(j).dgks_accept = accepted;
        diaginfo(j).final_norm = norm_new;
        diaginfo(j).norm_ratio = safe_ratio(norm_new,base);

        if ~(isfinite(norm_new) && norm_new > opts.orthTol*base)
            diaginfo(j).reason = 'relative_dependent_after_reorth';
            continue;
        end

        V(:,end+1) = z/norm_new; 
        added = added+1;
        diaginfo(j).appended = true;
        if accepted
            diaginfo(j).reason = 'accepted_relative_dgks';
        else
            diaginfo(j).reason = 'accepted_relative_after_max_reorth';
        end
    end
end

function r = safe_ratio(a,b)
if isfinite(a) && isfinite(b) && b > 0
    r = a/b;
else
    r = NaN;
end
end
