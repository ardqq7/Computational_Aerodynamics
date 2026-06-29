clear, clc,clf,close all;

% Mesh
[x, y] = meshgrid(-2:0.01:2, -2:0.01:2);
k = 1;

% Velocity
u = (k.*y) ./ (x.^2 + y.^2);
v = ( -k.*x ) ./ (x.^2 + y.^2);

speed = hypot(u,v);
contourf(x, y, speed, 20, 'LineStyle', 'none');
colormap(turbo);
colorbar;
hold on;

% Starting points for streamlines
startx = -5:0.1:5;
starty = -5:0.1:5;

% Streamlines
streamline(x, y, u, v, startx, starty);

xlabel('x', 'FontWeight', 'bold');
    ylabel('y', 'FontWeight', 'bold');
    title('Velocity Field', 'FontWeight', 'bold');
    axis equal, grid minor;
    xlim([-6,6]), ylim([-6,6]);
