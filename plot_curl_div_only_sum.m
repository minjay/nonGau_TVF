function plot_curl_div(theta, phi, theta_grid, phi_grid, vec_field_curl, vec_field_div, vec_field_sum,...
    mag_curl_grid, mag_div_grid, mag_sum_grid)
% plot curl-free, div-free and sum of them
% mag and quivers

fig = figure;
% find the max value
% max(mag_curl_grid) is a row vector
cmax = max([max(mag_curl_grid) max(mag_div_grid) max(mag_sum_grid)]);
% 0.62+0.06 = 0.68
subplot('position', [0.1 0.68 0.7 0.28])
plot_mag_quivers(theta, phi, vec_field_curl(:, 1), vec_field_curl(:, 2), theta_grid, phi_grid, mag_curl_grid)
caxis([0 cmax])
title('Curl-free')
% 0.28+0.06 = 0.34
% 0.34+0.28 = 0.62
subplot('position', [0.1 0.34 0.7 0.28])
plot_mag_quivers(theta, phi, vec_field_div(:, 1), vec_field_div(:, 2), theta_grid, phi_grid, mag_div_grid)
caxis([0 cmax])
title('Div-free')
% 0+0.28 = 0.28
subplot('position', [0.1 0 0.7 0.28])
plot_mag_quivers(theta, phi, vec_field_sum(:, 1), vec_field_sum(:, 2), theta_grid, phi_grid, mag_sum_grid)
caxis([0 cmax])
title('Sum')
h = colorbar;
set(h, 'Position', [0.85 0.02 0.05 0.96]);
% set size
% tall
set(fig, 'Position', [0, 0, 600, 800]);

end