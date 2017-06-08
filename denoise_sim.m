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

% we need to know/have a good guess of how many levels it has
% whether the noise level is relatively large or small
% if the noise level is large, we need to use fewer levels
% otherwise, we need to use more levels
A_norm_sub = [A_norm(:, 1:sum(Npix(1:2))) A_norm(:, M+1:M+sum(Npix(1:2)))];
% A^{+}Y
beta_est_sub = pinv(A_norm_sub)*Y;
beta_est = pinv(A_norm)*Y;
resid_sub = Y-A_norm_sub*beta_est_sub;
resid = Y-A_norm*beta_est;
tau_est = std(resid_sub)
std(resid)
% MAD/0.6745
median(abs(resid_sub-median(resid_sub)))/0.6745
median(abs(resid-median(resid)))/0.6745

factors = 0.05:0.05:1.5;
n_factor = length(factors);
lambdas_bpdn = factors*tau_est*sqrt(2*log(2*M));
Y_recover = cell(n_factor, 1);
betas = cell(n_factor, 1);
AIC = zeros(n_factor, 1);
AICc = zeros(n_factor, 1);
MSE = zeros(n_factor, 1);
df_zou = zeros(n_factor, 1);
% this is more appropriate
df_tib = zeros(n_factor, 1);
opts = as_setparms;
opts.loglevel = 1;
inform = [];
for i = 1:n_factor
    % "inform" is used for warm start
    [beta, inform] = as_bpdn(A_norm, Y, lambdas_bpdn(i), opts, inform);
    Y_recover{i} = A_norm*beta;
    betas{i} = beta;
    df_zou(i) = sum(beta~=0);
    df_tib(i) = rank(A_norm(:, beta~=0));
    AIC(i) = sum((Y-A_norm*beta).^2)/N/(tau_est^2)+2/N*df_tib(i);
    AICc(i) = AIC(i)+2*df_tib(i)*(df_tib(i)+1)/(N-df_tib(i)-1)/N;
    SSE(i) = sum((Y_recover{i}-Y_signal).^2);
end

% plot information criteria
figure
plot(factors, [SSE AIC AICc], '-o', 'LineWidth', 1.5)
legend('SSE', 'AIC', 'AICc', 'Location', 'Best')
set(gca, 'FontSize', 12)
hold on
[~, index_SSE] = min(SSE);
plot(factors(index_SSE), MSE(index_SSE), 'r*')
[~, index_AIC] = min(AIC);
plot(factors(index_AIC), AIC(index_AIC), 'r*')
xlabel('c_i')

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.05 0.05], [0.01 0.01]);
% between 3 and 5
scale = 4;
subplot(1, 2, 1)
plot_quivers_true_size(theta, phi, Y(1:N)/scale, Y(N+1:end)/scale, 'b')
hold on
draw_contour
title('Observed (Noisy)')
subplot(1, 2, 2)
plot_quivers_true_size(theta, phi, Y_signal(1:N)/scale, Y_signal(N+1:end)/scale, 'b')
hold on
plot_quivers_true_size(theta, phi, Y_recover{index_AIC}(1:N)/scale, Y_recover{index_AIC}(N+1:end)/scale, [1.000000 0.000000 0.500000])
draw_contour
title('Signal and Recovered')
legend('Signal', 'Recovered', 'Location', 'Best')
legend('boxoff')
