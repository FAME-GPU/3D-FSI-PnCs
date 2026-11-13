function Bi = Strain_matrix_tensor(a)
    
    b_val = a(:,1,:); 
    c_val = a(:,2,:); 
    d_val = a(:,3,:); 
    
    B_base = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0];
    C_base = [0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0];
    D_base = [0 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 0];
    
    Bi = b_val .* B_base + c_val .* C_base + d_val .* D_base;
end