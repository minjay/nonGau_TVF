function draw_contour

% draw contour
th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');
lam = pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');

end