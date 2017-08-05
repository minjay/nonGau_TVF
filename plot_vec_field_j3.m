clear
clc

B = 2;
j_min = 3;
j_max = 3;

% the sampling locations
Nside = 16;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end
res = 100;
% avoid the north and south poles
% phi ranges from -pi to pi
theta_grid = linspace(0+1e-3, pi-1e-3, res/2);
phi_grid = linspace(-pi, pi, res);
[phi_grid, theta_grid] = meshgrid(phi_grid, theta_grid);
theta_vec = theta_grid(:);
phi_vec = phi_grid(:);
N_grid = length(theta_vec);
% A_2_1 is the design matrix corresponding to the u component of Psi^{(2)}
% A_2_2 ... the v component of Psi^{(2)}
% A_3_1 ... the u component of Psi^{(3)}
% A_3_2 ... the v component of Psi^{(3)}
[Npix, grid_points, A_2_1, A_2_2, A_3_1, A_3_2] = get_A_curl_div_ss(B, j_min, j_max, theta, phi);
[~, ~, A_2_1_grid, A_2_2_grid, A_3_1_grid, A_3_2_grid] = get_A_curl_div_ss(B, j_min, j_max, theta_vec, phi_vec);

M = size(A_2_1, 2);
alpha = 3;
sigma = B.^((-alpha/2)*(j_min:j_max));
% highly correlated case
C = [1, 0.5; 0.5, 1];
nu = 2.5;
% set seed
rng('default')
gamma = get_coef_gamma(M, sigma, j_min, j_max, Npix, C, nu);
[vec_field_curl, vec_field_div, vec_field_sum,...
    mag_curl_grid, mag_div_grid, mag_sum_grid] = get_vec_field(N, N_grid, gamma,...
    A_2_1, A_2_2, A_3_1, A_3_2, A_2_1_grid, A_2_2_grid, A_3_1_grid, A_3_2_grid,...
    theta_grid);

% plot
plot_curl_div(theta, phi, theta_grid, phi_grid, vec_field_curl, vec_field_div, vec_field_sum,...
    mag_curl_grid, mag_div_grid, mag_sum_grid)
