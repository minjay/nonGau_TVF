function gamma = get_coef_gamma(M, sigma, j_min, j_max, Npix, C, nu)
% simulate coefficients gamma
% M: number of basis functions
% C: correlation matrix

gamma = zeros(M, 2);

st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    gamma(range, :) = mvtrnd(C, nu, Npix(index_j))*sigma(index_j);
    st = st+Npix(index_j);
end

end