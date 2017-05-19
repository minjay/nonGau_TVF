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
j_max = 3;
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

factors = [0.01 0.1 0.2 0.5];
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
scale = 100;
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
theta_tran = pi/2-theta;
phi_tran = phi;
phi_tran(phi>pi)=phi(phi>pi)-2*pi;

[HX, HY] = sph2hammer(phi_tran, theta_tran);
subplot(2, 1, 1)
scatter(HX, HY, [], sqrt(Y_signal(1:N).^2+Y_signal(N+1:end).^2), 'filled')
hold on
draw_contour
subplot(2, 1, 2)
scatter(HX, HY, [], sqrt(Y_recover{2}(1:N).^2+Y_recover{2}(N+1:end).^2), 'filled')
hold on
draw_contour