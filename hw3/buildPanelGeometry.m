function [geom, xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y)
% Build midpoints, tangents, normals from a closed contour (x,y).
% Returns fields: midpp (Nx2), S (Nx1), tx, ty, nx, ny, beta (Nx1), N

    % column vectors
    x = x(:); y = y(:);

    % Panel count
    N = numel(x) - 1; % -1 because start point and end point are same

    % Panel vectors and lengths
    dx = diff(x); 
    dy = diff(y);
    S  = hypot(dx, dy);

    % Unit tangents (from node i to i+1)
    tx = dx ./ S; 
    ty = dy ./ S;

    % Outward unit normals for CCW contour (rotate +90°)
    nx = -ty;
    ny =  tx;

    % Midpoints (collocation points)
    xc = 0.5*(x(1:end-1) + x(2:end));
    yc = 0.5*(y(1:end-1) + y(2:end));
    midpp = [xc, yc];

    % Panel angle wrt x-axis
    beta = atan2(dy, dx);


    % Normal angle (phi = beta + 90°)
    phi = beta + pi/2;     % radians

    % Pack results into struct
    geom = struct('N',N,'midpp',midpp,'S',S,'tx',tx,'ty',ty,...
                  'nx',nx,'ny',ny,'beta',beta,'phi',phi,'x',x,'y',y);

                      % ---- Display summary ----
    fprintf('\nPanel geometry:\n');
    fprintf('---------------------------------------------\n');
    fprintf('%5s %12s %12s %12s\n', 'i', 'beta (deg)', 'phi (deg)', 'Length');
    fprintf('---------------------------------------------\n');
    for i = 1:N
        fprintf('%5d %12.3f %12.3f %12.4f\n', ...
            i, rad2deg(beta(i)), rad2deg(phi(i)), S(i));
    end
    fprintf('---------------------------------------------\n');
end