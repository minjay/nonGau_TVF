function mag = mag_vec_need(B, J, dist)
% function that calculates the magnitude of vectorial needlets

N = length(dist);
mag = zeros(J+1, N);
l_min = ceil(B^(-1));
l_max = floor(B^(J+1));
P_l1 = cell(l_max-l_min+1, 1);
% precompute P_l1
for l = l_min:l_max
    tmp = legendre(l, dist);
    % m = 1
    P_l1{l-l_min+1} = tmp(2, :);
end
for j = 0:J
    j
    l_st = ceil(B^(j-1));
    l_en = floor(B^(j+1));
    for l = l_st:l_en
        mag(j+1, :) = mag(j+1, :)+fun_b(l/B^j, B)*(2*l+1)/4/pi*P_l1{l-l_min+1};
    end
    % lambda_j ~ B^{-2j}
    mag(j+1, :) = abs(mag(j+1, :)/(B^j));
end

end