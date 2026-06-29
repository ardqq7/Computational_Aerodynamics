function [sigma, A, b] = assembleSourceSystem(Uinf, alphaDeg, geom)
% Build and solve A*sigma = b for constant-source panel method
% BC: no-penetration at collocation points (midpoints)
% INPUTS:
%   Uinf, alphaDeg : free-stream speed and AoA (deg)
%   geom : struct with fields
%          x,y (closed), midpp(:,1/2)=xc,yc, S, tx,ty, nx,ny, beta, N
% OUTPUTS:
%   sigma : source strengths per panel (Nx1)
%   A, b  : influence matrix and RHS

    % Unpack geometry
    x    = geom.x;    y    = geom.y;
    xc   = geom.midpp(:,1);  yc = geom.midpp(:,2);
    S    = geom.S;
    nx   = geom.nx;   ny = geom.ny;
    beta = geom.beta;
    N    = geom.N;

    % Panel start nodes
    xb = x(1:end-1);  yb = y(1:end-1);

    % Free-stream components
    alpha = deg2rad(alphaDeg);
    Ufx = Uinf*cos(alpha); 
    Ufy = Uinf*sin(alpha);

    % Influence matrix
    A = zeros(N,N);
    b = zeros(N,1);

    % Build A from source-panel induced normal velocity
    % Using standard 2D constant-source panel formulas (Katz & Plotkin)
    eps0 = 1e-14;
    for i = 1:N
        for j = 1:N
            % Transform collocation i to panel-j local frame
            cb = cos(beta(j)); sb = sin(beta(j));
            xij = (xc(i)-xb(j))*cb + (yc(i)-yb(j))*sb;
            yij = -(xc(i)-xb(j))*sb + (yc(i)-yb(j))*cb;
            % avoid singular y
            if abs(yij) < eps0, yij = sign(yij + (yij==0))*eps0; end

            x1 = 0; x2 = S(j);
            % Influence of unit source on local velocities at (xij,yij)
            term1 = atan2((x2 - xij)*yij, (x2 - xij)^2 + yij^2) ...
                  - atan2((x1 - xij)*yij, (x1 - xij)^2 + yij^2);
            term2 = 0.5*log(((x2 - xij)^2 + yij^2)/((x1 - xij)^2 + yij^2));

            u_loc =  term1/(2*pi);
            v_loc =  term2/(2*pi);

            % Rotate back to global
            ug =  u_loc*cb - v_loc*sb;
            vg =  u_loc*sb + v_loc*cb;

            % Normal component at panel i
            A(i,j) = ug*nx(i) + vg*ny(i);
        end

        % RHS: - (U_free · n_i)
        b(i) = - (Ufx*nx(i) + Ufy*ny(i));
    end

    % Add closure condition: sum(sigma) = 0
    A(N,:) = 1.0;
    b(N)   = 0.0;

    % Solve
    sigma = A \ b;
end