function [n, area] = Outer_normal_vector_and_area_tensor(node_coords, v_ref)
    
    v12 = node_coords(:,2,:) - node_coords(:,1,:);
    v13 = node_coords(:,3,:) - node_coords(:,1,:);
    
    n_cross = cross(v12, v13, 1);
    area = sqrt(sum(n_cross.^2, 1)) / 2;
    n = n_cross ./ (2 * area);
    
    mask = sum(n .* v_ref, 1) < 0;
    n(:,:,mask) = -n(:,:,mask);
end