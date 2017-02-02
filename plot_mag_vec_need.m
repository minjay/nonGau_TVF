% script to plot mag_vec_need

% high-resolution
dist = 0:0.0001:1;
B = 2;
J = 4;

mag = mag_vec_need(B, J, dist);

% convert to angular/great-circle distance
angle = acos(dist);

plot([-angle fliplr(angle)], [mag fliplr(mag)], 'LineWidth', 1.5)
legend('j=0', 'j=1', 'j=2', 'j=3', 'j=4', 'Location', 'Best')
axis tight
xlabel('Great-circle distance')
ylabel('Magnitude of \Psi^{(i)}_{jk}, i=2,3')
set(gca, 'FontSize', 12)
