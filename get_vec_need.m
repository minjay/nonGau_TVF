function [Psi_2_uv, Psi_3_uv, theta, phi] = get_vec_need(j, B, Nside)
% Nside determines the size of the grid

% the sampling locations
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
% location of needlet is in the center
theta_need = pi/2;
phi_need = 0;
% convert to R^3
[x_need, y_need, z_need] = trans_coord(theta_need, phi_need);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end
[x, y, z] = trans_coord(theta, phi);
% compute inner products
inn_prod = x_need*x+y_need*y+z_need*z;
% check whether inner product==1 or -1
find(inn_prod==1)
find(inn_prod==-1)

l_st = ceil(B^(j-1));
l_en = floor(B^(j+1));
P_l1 = cell(l_en-l_st+1, 1);
% precompute P_l1
for l = l_st:l_en
    tmp = legendre(l, inn_prod);
    % m = 1
    P_l1{l-l_st+1} = tmp(2, :);
end

coef = zeros(1, N);
for l = l_st:l_en
    coef = coef-fun_b(l/B^j, B)*(2*l+1)/4/pi*P_l1{l-l_st+1}./sqrt(1-inn_prod'.^2);
end
vec_2 = zeros(3, N);
for i = 1:N
    vec_2(:, i) = [x_need; y_need; z_need]-inn_prod(i)*[x(i); y(i); z(i)];
end
vec_3 = zeros(3, N);
for i = 1:N
    vec_3(:, i) = cross([x_need; y_need; z_need], [x(i); y(i); z(i)]);
end

% lambda_j ~ B^{-2j}
Psi_2 = repmat(coef, 3, 1).*vec_2/B^j;
Psi_3 = repmat(coef, 3, 1).*vec_3/B^j;

% convert to (u, v) for plotting
Psi_2_uv = zeros(2, N);
Psi_3_uv = zeros(2, N);
for i = 1:N
    T = [-sin(phi(i)) cos(phi(i)) 0;...
        -cos(theta(i))*cos(phi(i)) -cos(theta(i))*sin(phi(i)) sin(theta(i))];
    Psi_2_uv(:, i) = T*Psi_2(:, i);
    Psi_3_uv(:, i) = T*Psi_3(:, i);
end

end