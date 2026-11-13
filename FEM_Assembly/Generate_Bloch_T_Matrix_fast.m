function [T_s, T_f] = Generate_Bloch_T_Matrix_fast(meshData, k_vector, A)

    np = size(meshData.nodes, 2);
    C  = A \ meshData.nodes;       
    tol = 1e-5;

    all_boundary = unique([meshData.Faces.left, meshData.Faces.right, ...
                           meshData.Faces.front, meshData.Faces.back, ...
                           meshData.Faces.bottom, meshData.Faces.top]);
    all_boundary = all_boundary(:);

    theta = [dot(k_vector, A(:,1));
             dot(k_vector, A(:,2));
             dot(k_vector, A(:,3))];

    is_slave = false(np, 1);
    cols = (1:np).';
    vals = ones(np, 1);

    Cb = C(:, all_boundary);
    shift_b = Cb > (0.5 - tol);       
    slave_mask_b = any(shift_b, 1);
    slave_nodes = all_boundary(slave_mask_b);

    if ~isempty(slave_nodes)
        is_slave(slave_nodes) = true;

        shift_s = double(shift_b(:, slave_mask_b));   
        cm = C(:, slave_nodes) - shift_s;             
        phase = theta.' * shift_s;                    

        keyTol = tol;
        Ckey  = int64(round(C.'  / keyTol));          
        cmKey = int64(round(cm.' / keyTol));          

        [key_unique, ia] = unique(Ckey, 'rows', 'stable');
        node_unique = ia(:);

        [tf, loc] = ismember(cmKey, key_unique, 'rows');
        master_idx = zeros(numel(slave_nodes), 1);
        master_idx(tf) = node_unique(loc(tf));

        miss = find(~tf);
        for jj = miss(:).'
            d2 = sum(bsxfun(@minus, C, cm(:, jj)).^2, 1);
            [d2min, master_idx(jj)] = min(d2);
            if d2min > (10 * tol)^2
                warning('Generate_Bloch_T_Matrix_fast:MasterFar', ...
                    'Slave node %d: nearest master distance = %.3e, larger than expected.', ...
                    slave_nodes(jj), sqrt(d2min));
            end
        end

        cols(slave_nodes) = master_idx;
        vals(slave_nodes) = exp(1i * phase(:));
    end

    f_nodes = meshData.nodeidx_f(:);
    f_reduced = f_nodes(~is_slave(f_nodes));

    [ok_f, col_f] = ismember(cols(f_nodes), f_reduced);
    row_f = (1:numel(f_nodes)).';
    if any(~ok_f)
        warning('Generate_Bloch_T_Matrix_fast:FluidMasterOutside', ...
            '%d fluid rows map to nodes outside f_reduced; these rows are left zero, same effect as submatrix extraction.', ...
            nnz(~ok_f));
    end
    T_f = sparse(row_f(ok_f), col_f(ok_f), vals(f_nodes(ok_f)), ...
                 numel(f_nodes), numel(f_reduced));

    s_nodes = meshData.nodeidx_s(:);
    s_reduced = s_nodes(~is_slave(s_nodes));

    [ok_s, col_s] = ismember(cols(s_nodes), s_reduced);
    row_s = (1:numel(s_nodes)).';
    if any(~ok_s)
        warning('Generate_Bloch_T_Matrix_fast:SolidMasterOutside', ...
            '%d solid rows map to nodes outside s_reduced; these rows are left zero, same effect as submatrix extraction.', ...
            nnz(~ok_s));
    end
    T_s_node = sparse(row_s(ok_s), col_s(ok_s), vals(s_nodes(ok_s)), ...
                      numel(s_nodes), numel(s_reduced));

    T_s = kron(T_s_node, speye(3));
end
