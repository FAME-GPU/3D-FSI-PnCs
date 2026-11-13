function [ phys ] = Physical_settings(model_type)

phys.model_type = model_type;

E = 210e9/1e10;   
nu = 0.3;   
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
mu = E / (2 * (1 + nu));
C_matrix = diag([lambda+2*mu, lambda+2*mu, lambda+2*mu, mu, mu, mu]);
C_matrix(1,2) = lambda; C_matrix(1,3) = lambda;
C_matrix(2,1) = lambda; C_matrix(2,3) = lambda;
C_matrix(3,1) = lambda; C_matrix(3,2) = lambda;
phys.meteria.c = C_matrix;
phys.meteria.rho_s = 7850 * eye(3); phys.meteria.rho_f_inv = (1 / 1000) * eye(3); 
phys.meteria.k_f = 1000 * 1500^2/1e10; 

switch phys.model_type
    case {'SC', 'SC_InnerFluid', 'SC_Quad'}
        phys.lattice.lattice_vector = [1 0 0;0 1 0;0 0 1];
end

B = 2*pi*inv(phys.lattice.lattice_vector');

switch phys.model_type
    case {'SC', 'SC_InnerFluid', 'SC_Quad'}
        path_string = 'GXMGRX|MR';
        vertex.G  = B*[    0,   0,   0 ]';
        vertex.M  = B*[  0.5,  0.5, 0.0 ]';
        vertex.R  = B*[  0.5,  0.5, 0.5 ]';
        vertex.X  = B*[    0,  0.5,   0 ]';
end

phys.lattice.path_string = path_string;
phys.lattice.vertex = vertex;

n = length(path_string);
part_num = 10;
phys.lattice.part_num = part_num;
wave_vec_array = [];
for i = 1:n-1
    if strcmp(path_string(i),'|') == 0 && strcmp(path_string(i+1),'|') == 0
        string_x = [ 'subpath(1,:) = linspace( vertex.',path_string(i),'(1), vertex.',path_string(i+1),'(1),', num2str(part_num),');'];
        string_y = [ 'subpath(2,:) = linspace( vertex.',path_string(i),'(2), vertex.',path_string(i+1),'(2),', num2str(part_num),');'];
        string_z = [ 'subpath(3,:) = linspace( vertex.',path_string(i),'(3), vertex.',path_string(i+1),'(3),', num2str(part_num),');'];
        eval(string_x); eval(string_y); eval(string_z);
        if i == n-1 || strcmp(path_string(i+2),'|') == 1
            wave_vec_array = [wave_vec_array,subpath(:,1:end)];
        else
            wave_vec_array = [wave_vec_array,subpath(:,1:end-1)];
        end
    end
end
phys.lattice.wave_vec_array = wave_vec_array;

end