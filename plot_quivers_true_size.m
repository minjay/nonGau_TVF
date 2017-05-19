function plot_quivers_true_size(theta, phi, u, v, color)
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
theta = pi/2-theta;
phi(phi>pi)=phi(phi>pi)-2*pi;

[HX, HY] = sph2hammer(phi, theta);
quiver(HX, HY, u, v, 'color', color, 'Autoscale', 'off');

axis equal
axis tight
axis off

end