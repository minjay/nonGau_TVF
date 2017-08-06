clear
clc

subplot = @(m,n,p) subtightplot (m, n, p, [0 0.05], [0 0.02], [0.05 0.01]);

j = 2;
B = 2;
Nside = 16;
[Psi_2_uv_j2, Psi_3_uv_j2, theta, phi] = get_vec_need(j, B, Nside);
subplot(2, 2, 1)
plot_quivers(theta, phi, Psi_2_uv_j2(1, :)', Psi_2_uv_j2(2, :)')
title('i=2, j=2')
subplot(2, 2, 2)
plot_quivers(theta, phi, Psi_3_uv_j2(1, :)', Psi_3_uv_j2(2, :)')
title('i=3, j=2')

j = 3;
[Psi_2_uv_j3, Psi_3_uv_j3, theta, phi] = get_vec_need(j, B, Nside);
subplot(2, 2, 3)
plot_quivers(theta, phi, Psi_2_uv_j3(1, :)', Psi_2_uv_j3(2, :)')
title('i=2, j=3')
subplot(2, 2, 4)
plot_quivers(theta, phi, Psi_3_uv_j3(1, :)', Psi_3_uv_j3(2, :)')
title('i=3, j=3')

theta = pi/2-theta;
phi(phi>pi)=phi(phi>pi)-2*pi;

[HX, HY] = sph2hammer(phi, theta);

th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
th2 = th;
lam2 = pi+0*th;
[xh2,yh2] = sph2hammer(lam2,th2);

save('vec_need.mat', 'HX', 'HY', 'th', 'lam', 'th2', 'lam2', 'Psi_2_uv_j2', 'Psi_3_uv_j2',...
    'Psi_2_uv_j3', 'Psi_3_uv_j3')