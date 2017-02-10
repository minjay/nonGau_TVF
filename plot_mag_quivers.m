function plot_mag_quivers(theta, phi, u, v, theta_grid, phi_grid, f)
%PLOT_QUIVERS   Plots simulated samples on the sphere (after Hammer
%   projection).
%
%   PLOT_QUIVERS(THETA, PHI, U, V);
%
% Inputs:
%   THETA, PHI - the sampling locations
%   U, V - the u and v components
%
% Author: Minjie Fan, 2015

% convert theta and phi
theta_grid = pi/2-theta_grid;
theta = pi/2-theta;
phi(phi>pi) = phi(phi>pi)-2*pi;
[HX_grid, HY_grid] = sph2hammer(phi_grid, theta_grid);
[HX, HY] = sph2hammer(phi, theta);
pcolor(HX_grid, HY_grid, f);
shading flat
hold on
quiver(HX, HY, u, v, 'k');

% draw contour
th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');
lam = pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');

axis equal
axis tight
axis off

set(gca, 'FontSize', 12)

end