function [x, y, xc, yc, tx, ty, nx, ny, S] = makeCirclePanels(r, N, cx, cy, startAngleDeg, fileName)

    % Angles (CCW) for nodes; ensure closure by repeating first node at end
    theta0 = deg2rad(startAngleDeg);
    dth    = 2*pi/N;
    theta  = theta0 + (0:N)*dth;   % length N+1

    % Node coordinates (closed contour)
    x = cx + r*cos(theta).';
    y = cy + r*sin(theta).';

    % Panel vectors and lengths
    dx = diff(x);                   % length N
    dy = diff(y);
    S  = hypot(dx,dy);

    % Unit tangents (CCW)
    tx = dx./S;
    ty = dy./S;

    % Outward unit normals (rotate tangent +90°)
    nx = -ty;
    ny =  tx;

    % Collocation (midpoint) coordinates
    xc = 0.5*(x(1:end-1) + x(2:end));
    yc = 0.5*(y(1:end-1) + y(2:end));

    % Optional: write to CSV
    if strlength(fileName) > 0
        T = table(x, y);
        writetable(T, fileName);
    end
end