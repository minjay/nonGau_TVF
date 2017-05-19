clear
clc

addpath(genpath('/Users/minjay/Documents/MATLAB/asp'))
% set seed
rng(1)

B = 2;
j_min = 1;
j_max = 3;

% the sampling locations
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end
% A_2_1 is the design matrix corresponding to the u component of Psi^{(2)}
% A_2_2 ... the v component of Psi^{(2)}
% A_3_1 ... the u component of Psi^{(3)}
% A_3_2 ... the v component of Psi^{(3)}
[Npix, grid_points, A_2_1, A_2_2, A_3_1, A_3_2] = get_A_curl_div_ss(B, j_min, j_max, theta, phi);

M = size(A_2_1, 2);
alpha = 3;
sigma = B.^((-alpha/2)*(j_min:j_max));
C = [1, 0; 0, 1];
nu = 3;
gamma = get_coef_gamma(M, sigma, j_min, j_max, Npix, C, nu);

A_all = [A_2_1 A_3_1; A_2_2 A_3_2];
% normalize cols of A
A_norm = zeros(size(A_all));
for i = 1:M*2
    A_norm(:, i) = A_all(:, i)./norm(A_all(:, i));
end

gamma_all = [gamma(:, 1); gamma(:, 2)];

% white noise std
% between 0.1 and 0.2
tau = 0.15;
Y_signal = A_norm*gamma_all;
Y = Y_signal+tau*randn(N*2, 1);

factors = [0.1 0.2 0.5 1];
n_factor = length(factors);
lambdas_bpdn = factors*tau*sqrt(2*log(2*M));
Y_recover = cell(n_factor, 1);
betas = cell(n_factor, 1);
for i = 1:n_factor
    opts = as_setparms;
    opts.loglevel = 1;
    inform = [];
    [beta, inform] = as_bpdn(A_norm, Y, lambdas_bpdn(i), opts, inform);
    Y_recover{i} = A_norm*beta;
    betas{i} = beta;
end

figure
% between 3 and 5
scale = 4;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.05 0.05], [0.01 0.01]);

subplot(3, 2, 1)
plot_quivers_true_size(theta, phi, Y(1:N)/scale, Y(N+1:end)/scale, 'b')
hold on
draw_contour
title('Observed (Noisy)')

subplot(3, 2, 2)
plot_quivers_true_size(theta, phi, Y_signal(1:N)/scale, Y_signal(N+1:end)/scale, 'b')
hold on
draw_contour
title('Signal')

for i = 1:n_factor
    subplot(3, 2, i+2)
    plot_quivers_true_size(theta, phi, Y_recover{i}(1:N)/scale, Y_recover{i}(N+1:end)/scale, 'b')
    hold on
    draw_contour
    title(['Recovered (alpha=', num2str(factors(i)), ')'])
end

figure
plot_quivers_true_size(theta, phi, Y_signal(1:N)/scale, Y_signal(N+1:end)/scale, 'b')
hold on
plot_quivers_true_size(theta, phi, Y_recover{2}(1:N)/scale, Y_recover{2}(N+1:end)/scale, [1.000000 0.000000 0.500000])
draw_contour
title('Signal and Recovered')
legend('Signal', 'Recovered', 'Location', 'bestoutside')
legend('boxoff')

figure
clear subplot
subplot(1, 2, 1)
stem(gamma(:, 1))
hold on
plot(betas{2}(1:M), 'r*')
legend('True', 'Recovered')
title('Curl-free')
subplot(1, 2, 2)
stem(gamma(:, 2))
hold on
plot(betas{2}(M+1:end), 'r*')
legend('True', 'Recovered')
title('Div-free')