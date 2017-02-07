clear
clc

subplot = @(m,n,p) subtightplot (m, n, p, [0 0.05], [0 0.02], [0.05 0.01]);

j = 2;
B = 2;
Nside = 16;
[Psi_2_uv, Psi_3_uv, theta, phi] = get_vec_need(j, B, Nside);
subplot(2, 2, 1)
plot_quivers(theta, phi, Psi_2_uv(1, :)', Psi_2_uv(2, :)')
title('i=2, j=2')
subplot(2, 2, 2)
plot_quivers(theta, phi, Psi_3_uv(1, :)', Psi_3_uv(2, :)')
title('i=3, j=2')

j = 3;
[Psi_2_uv, Psi_3_uv, theta, phi] = get_vec_need(j, B, Nside);
subplot(2, 2, 3)
plot_quivers(theta, phi, Psi_2_uv(1, :)', Psi_2_uv(2, :)')
title('i=2, j=3')
subplot(2, 2, 4)
plot_quivers(theta, phi, Psi_3_uv(1, :)', Psi_3_uv(2, :)')
title('i=3, j=3')
