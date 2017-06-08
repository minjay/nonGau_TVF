clear

addpath(genpath('/Users/minjay/Documents/MATLAB/asp'))

% read data
uwnd_all = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'uas');
% lat and lon (in radian)
% lat in [-pi/2, pi/2]
% lon in [0, 2*pi]
lat = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'lat')/180*pi;
lon = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'lon')/180*pi;
vwnd_all = ncread('vas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'vas');

% extract the data from Jan. 2000 to Dec. 2008
uwnd = uwnd_all(:, :, 6);
vwnd = vwnd_all(:, :, 6);

res = 5;
lat_sub = lat(1:res:end);
lon_sub = lon(1:res:end);
uwnd_sub = uwnd(1:res:end, 1:res:end);
vwnd_sub = vwnd(1:res:end, 1:res:end);

[lat_m, lon_m] = meshgrid(lat_sub, lon_sub);

% reshape
uwnd_sub = uwnd_sub(:);
vwnd_sub = vwnd_sub(:);
lat_m = lat_m(:);
lon_m = lon_m(:);

% find the indices of NaN
idx_NaN = isnan(uwnd_sub)+isnan(vwnd_sub);
uwnd_sub_v = uwnd_sub(~idx_NaN);
vwnd_sub_v = vwnd_sub(~idx_NaN);
lat_m_v = lat_m(~idx_NaN);
lon_m_v = lon_m(~idx_NaN);

theta = pi/2-lat_m_v;
phi = lon_m_v;

B = 2;
j_min = 0;
j_max = 4;
[Npix, grid_points, A_2_1, A_2_2, A_3_1, A_3_2] = get_A_curl_div_ss(B, j_min, j_max, theta, phi);
M = size(A_2_1, 2);

A_all = [A_2_1 A_3_1; A_2_2 A_3_2];
% normalize cols of A
A_norm = zeros(size(A_all));
for i = 1:M*2
    A_norm(:, i) = A_all(:, i)./norm(A_all(:, i));
end

Y_signal = [uwnd_sub_v; vwnd_sub_v];
tau = 2;
N = length(uwnd_sub_v);
% set seed
rng(1)
Y = Y_signal+tau*randn(N*2, 1);

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
L2_err = zeros(n_factor, 1);
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
    L2_err(i) = norm(Y_recover{i}-Y_signal);
end

% plot information criteria
figure
plot(factors, [L2_err AIC AICc], '-o', 'LineWidth', 1.5)
legend('L2-err', 'AIC', 'AICc', 'Location', 'Best')
set(gca, 'FontSize', 12)
hold on
[~, index_L2_err] = min(L2_err);
plot(factors(index_L2_err), L2_err(index_L2_err), 'r*')
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
