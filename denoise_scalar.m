% tests for scalar needlets

addpath(genpath('/Users/minjay/Documents/MATLAB/asp'))

clear
clc

rng(1)

nu = 3;
alpha = 3;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

j_min = 1;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
M = size(A, 2);

sigma_j = B.^(-alpha/2*(j_min:j_max));
sigma_j = sigma_j/sigma_j(1);

c = zeros(M, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    st = st+Npix(index_j);
end

% normalize cols of A
A_norm = zeros(size(A));
for i = 1:M
    A_norm(:, i) = A(:, i)./norm(A(:, i));
end

% white noise std
tau = 0.5;
Y_signal = A_norm*c;
Y = Y_signal+randn(N, 1)*tau;

factor = 0.2;
lambda_bpdn = factor*tau*sqrt(2*log(M));
opts          = as_setparms;
opts.loglevel = 1;
inform        = [];  % IMPORTANT: must initialize in this way.
[beta, inform] = as_bpdn(A_norm, Y, lambda_bpdn, opts, inform);
Y_recover = A_norm*beta;

h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.02], [0.1 0.1], [0.075 0.01]);
cmin = min([min(Y) min(Y_signal) min(Y_recover)]);
cmax = max([max(Y) max(Y_signal) max(Y_recover)]);
subplot(1, 3, 1)
scatter(theta, phi, [], Y_signal, 'filled')
axis tight
colorbar
caxis([cmin cmax])
title('Signal')
subplot(1, 3, 2)
scatter(theta, phi, [], Y, 'filled')
axis tight
colorbar
caxis([cmin cmax])
title('Noisy')
subplot(1, 3, 3)
scatter(theta, phi, [], Y_recover, 'filled')
axis tight
colorbar
caxis([cmin cmax])
title('Recovered')
set(h, 'Position', [0, 0, 500, 200]);

figure
scatter(Y_signal, Y_recover)
axis equal
axis tight
hline = refline(1, 0);
hline.Color='r';
xlabel('Signal')
ylabel('Recovered')

norm(Y_signal-Y_recover)

figure
stem(c)
hold on
plot(beta, 'r*')
y_range = get(gca, 'ylim');
x_right = cumsum(Npix)+0.5;
for i = 1:(length(x_right)-1)
    plot([x_right(i) x_right(i)], y_range, '--', 'Color', [0.6 0.6 0.6])
end
legend('True', 'Recovered')
axis tight

% check sparsity of beta
figure
clear subplot
subplot(1, 2, 1)
histogram(c, 100)
title('True')
subplot(1, 2, 2)
histogram(beta, 100)
title('Recovered')