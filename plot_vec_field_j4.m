clear
clc

B = 2;
j_min = 4;
j_max = 4;

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

% sanity checking
[gamma_t, gamma_g] = get_coef_gamma(1e7, sigma, j_min, j_max, 1e7, C, nu);
% std(gamma_t) is more unstable
std(gamma_t)
std(gamma_g)
sigma*sqrt(nu/(nu-2))

% set seed
rng(1)
[gamma_t, gamma_g] = get_coef_gamma(M, sigma, j_min, j_max, Npix, C, nu);


[vec_field_curl_t, vec_field_div_t, vec_field_sum_t,...
    mag_curl_grid_t, mag_div_grid_t, mag_sum_grid_t] = get_vec_field(N, N_grid, gamma_t,...
    A_2_1, A_2_2, A_3_1, A_3_2, A_2_1_grid, A_2_2_grid, A_3_1_grid, A_3_2_grid,...
    theta_grid);

[vec_field_curl_g, vec_field_div_g, vec_field_sum_g,...
    mag_curl_grid_g, mag_div_grid_g, mag_sum_grid_g] = get_vec_field(N, N_grid, gamma_g,...
    A_2_1, A_2_2, A_3_1, A_3_2, A_2_1_grid, A_2_2_grid, A_3_1_grid, A_3_2_grid,...
    theta_grid);

theta = pi/2-theta;
phi(phi>pi) = phi(phi>pi)-2*pi;
[HX, HY] = sph2hammer(phi, theta);

th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
th2 = th;
lam2 = pi+0*th;
[xh2,yh2] = sph2hammer(lam2,th2);

save('vec_field_j4.mat', 'HX', 'HY', 'vec_field_sum_t', 'vec_field_sum_g',...
    'xh', 'yh', 'xh2', 'yh2')