function matrices = Assembly_mass_and_stiffness_matrix(meshData, meterial)

    coords   = meshData.nodes;             
    elements = meshData.elements;        
    matID    = meshData.cellIndices;
    num_node = size(coords, 2);                
    
    template3 = (1/12) * [2 1 1; 1 2 1; 1 1 2];
    template4 = (1/20) * [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2];
    
    rho_s     = meterial.rho_s;      
    rho_f_inv = meterial.rho_f_inv;  
    k_f       = meterial.k_f;        
    c_tensor  = meterial.c;          
    
    elem_idx_s = find(matID == meshData.solidTag);
    elem_idx_f = find(matID == meshData.fluidTag);
    num_f = length(elem_idx_f);
    num_s = length(elem_idx_s);

    node_idx_f  = elements(:, elem_idx_f);
    node_coor_f = reshape(coords(:, node_idx_f), 3, 4, []);
    
    v12_f = node_coor_f(:,2,:) - node_coor_f(:,1,:);
    v13_f = node_coor_f(:,3,:) - node_coor_f(:,1,:);
    v14_f = node_coor_f(:,4,:) - node_coor_f(:,1,:);
    vol_f = abs(sum(v12_f .* cross(v13_f, v14_f, 1), 1)) / 6; 
    
    P_inv_f = pageinv([ones(1, 4, num_f); node_coor_f]);         
    B_all_f = P_inv_f(:, 2:4, :);
    
    [idx_i_f, idx_j_f] = meshgrid(1:4, 1:4);
    I_f = reshape(node_idx_f(idx_i_f(:), :), [], 1);
    J_f = reshape(node_idx_f(idx_j_f(:), :), [], 1);
    
    M_f_Tensor = vol_f .* (1/k_f * template4);
    K_f_Tensor = vol_f .* pagemtimes(pagemtimes(B_all_f, rho_f_inv), pagetranspose(B_all_f)); 

    node_idx_s  = elements(:, elem_idx_s);
    node_coor_s = reshape(coords(:, node_idx_s), 3, 4, []);
    
    node_idx_s_3dof = [3*node_idx_s(1,:)-2; 3*node_idx_s(1,:)-1; 3*node_idx_s(1,:); ...
                       3*node_idx_s(2,:)-2; 3*node_idx_s(2,:)-1; 3*node_idx_s(2,:); ...
                       3*node_idx_s(3,:)-2; 3*node_idx_s(3,:)-1; 3*node_idx_s(3,:); ...
                       3*node_idx_s(4,:)-2; 3*node_idx_s(4,:)-1; 3*node_idx_s(4,:)];
                   
    v12_s = node_coor_s(:,2,:) - node_coor_s(:,1,:);
    v13_s = node_coor_s(:,3,:) - node_coor_s(:,1,:);
    v14_s = node_coor_s(:,4,:) - node_coor_s(:,1,:);
    vol_s = abs(sum(v12_s .* cross(v13_s, v14_s, 1), 1)) / 6;
    
    P_inv_s  = pageinv([ones(1, 4, num_s); node_coor_s]);         
    B_all_s  = P_inv_s(:, 2:4, :);
    B_Strain = [Strain_matrix_tensor(B_all_s(1,:,:)), Strain_matrix_tensor(B_all_s(2,:,:)), ...
                Strain_matrix_tensor(B_all_s(3,:,:)), Strain_matrix_tensor(B_all_s(4,:,:))]; 
            
    [idx_i_s, idx_j_s] = meshgrid(1:12, 1:12);
    I_s = reshape(node_idx_s_3dof(idx_i_s(:), :), [], 1);
    J_s = reshape(node_idx_s_3dof(idx_j_s(:), :), [], 1);
    
    M_s_Tensor = vol_s .* kron(template4, rho_s);
    K_s_Tensor = vol_s .* pagemtimes(pagetranspose(B_Strain), pagemtimes(c_tensor, B_Strain));

    node_idx_int = meshData.interfaceTriangles; 
    num_int      = size(node_idx_int, 2);
    node_coor_int = reshape(coords(:, node_idx_int), 3, 3, []);
    
    t_solid      = elements(:, elem_idx_s);   
    all_faces    = [t_solid([1 2 3],:), t_solid([1 2 4],:), t_solid([1 3 4],:), t_solid([2 3 4],:)];
    nodes4       = [t_solid(4,:),       t_solid(3,:),       t_solid(2,:),       t_solid(1,:)];
    
    [~, loc]     = ismember(sort(node_idx_int, 1)', sort(all_faces, 1)', 'rows');
    v_ref_idx    = nodes4(loc);
    
    v_ref_vector = reshape(coords(:, v_ref_idx), 3, 1, []) - node_coor_int(:,1,:);
    [n_vec, area_int] = Outer_normal_vector_and_area_tensor(node_coor_int, v_ref_vector);
    
    n_vec_t = pagetranspose(n_vec);
    C_Tensor = zeros(3, 9, num_int);
    for row = 1:3
        for col = 1:3
            C_Tensor(row, (col-1)*3+1 : col*3, :) = template3(row,col) .* n_vec_t;
        end
    end
    C_Tensor = reshape(area_int, 1, 1, num_int) .* C_Tensor;
    
    node_idx_int_s = [3*node_idx_int(1,:)-2; 3*node_idx_int(1,:)-1; 3*node_idx_int(1,:); ...
                      3*node_idx_int(2,:)-2; 3*node_idx_int(2,:)-1; 3*node_idx_int(2,:); ...
                      3*node_idx_int(3,:)-2; 3*node_idx_int(3,:)-1; 3*node_idx_int(3,:)];
                  
    [idx_i_int, idx_j_int] = ndgrid(1:3, 1:9);
    I_int = reshape(node_idx_int(idx_i_int(:), :), [], 1);
    J_int = reshape(node_idx_int_s(idx_j_int(:), :), [], 1);

    M_f = sparse(I_f, J_f, M_f_Tensor(:), num_node, num_node);
    K_f = sparse(I_f, J_f, K_f_Tensor(:), num_node, num_node);
    
    M_s = sparse(I_s, J_s, M_s_Tensor(:), 3*num_node, 3*num_node);
    K_s = sparse(I_s, J_s, K_s_Tensor(:), 3*num_node, 3*num_node);
    
    C_mat = sparse(I_int, J_int, C_Tensor(:), num_node, 3*num_node); 

    dofs_f = meshData.nodeidx_f(:); 
    idx_s  = meshData.nodeidx_s(:)'; 
    dofs_s = reshape([3*idx_s - 2; 3*idx_s - 1; 3*idx_s], [], 1); 
    
    matrices = struct(...
        'M_f_tight', M_f(dofs_f, dofs_f), ...
        'K_f_tight', K_f(dofs_f, dofs_f), ...
        'M_s_tight', M_s(dofs_s, dofs_s), ...
        'K_s_tight', K_s(dofs_s, dofs_s), ...
        'C_tight',   C_mat(dofs_f, dofs_s) ...
    );

end