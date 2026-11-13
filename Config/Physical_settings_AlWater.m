function phys = Physical_settings_AlWater(model_type)

if nargin < 1 || isempty(model_type)
    model_type = 'P4mmm';
end
phys.model_type = model_type;

phys.geometry.a0 = 3.0e-3;
phys.geometry.d0 = 2.0e-3;
phys.geometry.tm = 2.0e-3;
phys.geometry.tw = 2.0e-3;
phys.geometry.az = phys.geometry.tm + phys.geometry.tw;
phys.geometry.cell_size = [phys.geometry.a0, ...
                           phys.geometry.a0, ...
                           phys.geometry.az];

stiffness_scale = 1.0e10;

E   = 69.0e9 / stiffness_scale;
nu  = 0.33;
rho_s = 2700;

lambda = E * nu / ((1 + nu) * (1 - 2*nu));
mu     = E / (2*(1 + nu));

C_matrix = diag([lambda + 2*mu, lambda + 2*mu, lambda + 2*mu, ...
                 mu, mu, mu]);
C_matrix(1,2) = lambda; C_matrix(1,3) = lambda;
C_matrix(2,1) = lambda; C_matrix(2,3) = lambda;
C_matrix(3,1) = lambda; C_matrix(3,2) = lambda;

rho_w = 1000;
c_w   = 1490;
bulk_w = rho_w*c_w^2 / stiffness_scale;

phys.meteria.c         = C_matrix;
phys.meteria.rho_s     = rho_s * eye(3);
phys.meteria.rho_f_inv = (1/rho_w) * eye(3);
phys.meteria.k_f       = bulk_w;

phys.material.Aluminum.E_Pa   = 69.0e9;
phys.material.Aluminum.nu     = nu;
phys.material.Aluminum.rho    = rho_s;
phys.material.Water.rho       = rho_w;
phys.material.Water.c         = c_w;
phys.material.Water.bulk_Pa   = rho_w*c_w^2;
phys.scaling.stiffness        = stiffness_scale;

a0 = phys.geometry.a0;
az = phys.geometry.az;

phys.lattice.lattice_vector = [a0, 0,  0; ...
                               0,  a0, 0; ...
                               0,  0,  az];

B = 2*pi*inv(phys.lattice.lattice_vector');

vertex.G = B*[0.0; 0.0; 0.0];
vertex.X = B*[0.0; 0.5; 0.0];
vertex.M = B*[0.5; 0.5; 0.0];
vertex.Z = B*[0.0; 0.0; 0.5];
vertex.R = B*[0.5; 0.5; 0.5]; 
vertex.A = B*[0.0; 0.5; 0.5]; 

path_string = 'XGMXAZRAG';
part_num = 10;

phys.lattice.path_string = path_string;
phys.lattice.vertex = vertex;
phys.lattice.part_num = part_num;

wave_vec_array = [];
for i = 1:length(path_string)-1
    p0 = vertex.(path_string(i));
    p1 = vertex.(path_string(i+1));
    subpath = [linspace(p0(1), p1(1), part_num); ...
               linspace(p0(2), p1(2), part_num); ...
               linspace(p0(3), p1(3), part_num)];
    if i < length(path_string)-1
        wave_vec_array = [wave_vec_array, subpath(:,1:end-1)]; 
    else
        wave_vec_array = [wave_vec_array, subpath]; 
    end
end
phys.lattice.wave_vec_array = wave_vec_array;

phys.frequency.lambda_to_Hz = @(lambda) ...
    sqrt(max(real(lambda),0)*stiffness_scale)/(2*pi);
phys.frequency.lambda_to_MHz = @(lambda) ...
    sqrt(max(real(lambda),0)*stiffness_scale)/(2*pi*1e6);

end
