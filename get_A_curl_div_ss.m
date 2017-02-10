function [Npix, grid_points, A_2_1, A_2_2, A_3_1, A_3_2] = get_A_curl_div_ss(B, j_min, j_max, theta, phi)
% get design matrices with size # of sampling locations * # of needlets
% A_2_1 is the design matrix corresponding to the u component of Psi^{(2)}
% A_2_2 ... the v component of Psi^{(2)}
% A_3_1 ... the u component of Psi^{(3)}
% A_3_2 ... the v component of Psi^{(3)}

load('ss.mat')

N = length(theta);
[x, y, z] = trans_coord(theta, phi);
len_j = j_max-j_min+1;

Npix = zeros(len_j, 1);
grid_points = cell(len_j, 1);
for j = j_min:j_max
    index_j = j-j_min+1;
    % t needs to be odd
    % the quadrature formula is exact for all polynomials of degree<=t
    t = 2*floor(B^(j+1))+1;
    grid_points{index_j} = ss{degree_t==t};
    Npix(index_j) = size(grid_points{index_j}, 1);
end

M = sum(Npix);
A_2_1 = zeros(N, M);
A_2_2 = zeros(N, M);
A_3_1 = zeros(N, M);
A_3_2 = zeros(N, M);
index_col = 0;

for j = j_min:j_max
    disp(['j = ', num2str(j), ' starts...'])
    index_j = j-j_min+1;
    sqrt_lambda = sqrt(4*pi/Npix(index_j));
    for k = 1:Npix(index_j)
        index_col = index_col+1;
        xyz_xi = grid_points{index_j}(k, :);
        [Psi_2_uv, Psi_3_uv] = eval_vec_need(B, j, xyz_xi, x, y, z, theta, phi, sqrt_lambda);
        A_2_1(:, index_col) = Psi_2_uv(1, :);
        A_2_2(:, index_col) = Psi_2_uv(2, :);
        A_3_1(:, index_col) = Psi_3_uv(1, :);
        A_3_2(:, index_col) = Psi_3_uv(2, :);
    end
end

end