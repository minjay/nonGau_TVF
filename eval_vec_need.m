function [Psi_2_uv, Psi_3_uv] = eval_vec_need(B, j, xyz_xi, x, y, z, theta, phi, sqrt_lambda)
% evaluate vectorial needlets with centroid xyz_xi at location (x, y, z) or
% (theta, phi)
% results are returned in terms of u and v components

% compute inner products between xyz_xi and (x, y, z)
N = length(x);
dist = zeros(N, 1);
for i = 1:N
    dist(i) = sum([x(i) y(i) z(i)].*xyz_xi);
    dist(i) = min(dist(i), 1);
    dist(i) = max(dist(i), -1);
end

l_st = ceil(B^(j-1));
l_en = floor(B^(j+1));
P_l1 = cell(l_en-l_st+1, 1);
% precompute P_l1
for l = l_st:l_en
    tmp = legendre(l, dist);
    % m = 1
    P_l1{l-l_st+1} = tmp(2, :);
end

coef = zeros(1, N);
for l = l_st:l_en
    coef = coef-fun_b(l/B^j, B)*(2*l+1)/4/pi*P_l1{l-l_st+1}./sqrt(1-dist'.^2);
end
vec_2 = zeros(3, N);
for i = 1:N
    vec_2(:, i) = xyz_xi'-dist(i)*[x(i); y(i); z(i)];
end
vec_3 = zeros(3, N);
for i = 1:N
    vec_3(:, i) = cross(xyz_xi', [x(i); y(i); z(i)]);
end

Psi_2 = repmat(coef, 3, 1).*vec_2*sqrt_lambda;
Psi_3 = repmat(coef, 3, 1).*vec_3*sqrt_lambda;

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