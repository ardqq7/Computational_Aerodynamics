function plotPanelsMidTangNorm(geom, scaleTN)
% Plot panels, midpoints (midpp), and optional tangent/normal vectors.
% scaleTN: arrow scale for quiver (default 0.25)

    if nargin<2, scaleTN = 0.25; end

    x = geom.x; y = geom.y;
    xc = geom.midpp(:,1); yc = geom.midpp(:,2);

    figure(2);

    set(gcf, 'Color','w'); hold on; axis equal; grid on; box on;

    % Panels (black)
    plot(x, y, 'k-', 'LineWidth', 1.5);

    % Nodes (blue, excluding the repeated last)
    plot(x(1:end-1), y(1:end-1), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 5);

    % Midpoints (red)
    plot(xc, yc, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 6);

    % Tangents (gray) & normals (magenta)
    quiver(xc, yc, geom.tx, geom.ty, scaleTN, 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    quiver(xc, yc, geom.nx, geom.ny, scaleTN, 'm', 'LineWidth', 1);

    title(sprintf('Panels, midpp, tangents & normals (N=%d)', geom.N));
    xlabel('x'); ylabel('y');
    legend({'Panels','Nodes','midpp','tangent t','normal n'}, 'Location','best');
end