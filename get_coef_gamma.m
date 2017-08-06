function [gamma_t, gamma_g] = get_coef_gamma(M, sigma, j_min, j_max, Npix, C, nu)
% simulate coefficients gamma
% M: number of basis functions
% C: correlation matrix

gamma_t = zeros(M, 2);
gamma_g = zeros(M, 2);

st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    gamma_g(range, :) = mvnrnd(zeros(1, size(C, 1)), C, Npix(index_j))*sigma(index_j);
    denom = sqrt(chi2rnd(nu, Npix(index_j), 1)/nu);
    gamma_t(range, 1) = gamma_g(range, 1)./denom;
    gamma_t(range, 2) = gamma_g(range, 2)./denom;
    gamma_g(range, :) = gamma_g(range, :)*sqrt(nu/(nu-2));
    st = st+Npix(index_j);
end

end