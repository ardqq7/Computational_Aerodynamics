function [Ut, Cp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma)
% Compute tangential velocity at midpoints and Cp distribution
% INPUTS:
%   Uinf, alphaDeg : free-stream
%   geom           : geometry struct (same as above)
%   sigma          : source strengths
% OUTPUTS:
%   Ut : tangential velocity at each panel midpoint
%   Cp : surface pressure coefficient = 1 - (Ut/Uinf)^2

    x    = geom.x;    y    = geom.y;
    xb   = x(1:end-1); yb = y(1:end-1);
    S    = geom.S;    beta = geom.beta;
    tx   = geom.tx;   ty  = geom.ty;
    xc   = geom.midpp(:,1);  yc = geom.midpp(:,2);
    N    = geom.N;

    alpha = deg2rad(alphaDeg);
    Ufx = Uinf*cos(alpha); 
    Ufy = Uinf*sin(alpha);

    Ut = zeros(N,1);
    eps0 = 1e-14;

    for i = 1:N
        % Free-stream tangential at panel i
        U_t = Ufx*tx(i) + Ufy*ty(i);

        % Induced tangential by all source panels
        Uind_t = 0.0;
        for j = 1:N
            cb = cos(beta(j)); sb = sin(beta(j));
            xij = (xc(i)-xb(j))*cb + (yc(i)-yb(j))*sb;
            yij = -(xc(i)-xb(j))*sb + (yc(i)-yb(j))*cb;
            if abs(yij) < eps0, yij = sign(yij + (yij==0))*eps0; end

            x1 = 0; x2 = S(j);
            term1 = atan2((x2 - xij)*yij, (x2 - xij)^2 + yij^2) ...
                  - atan2((x1 - xij)*yij, (x1 - xij)^2 + yij^2);
            term2 = 0.5*log(((x2 - xij)^2 + yij^2)/((x1 - xij)^2 + yij^2));

            u_loc =  term1/(2*pi);
            v_loc =  term2/(2*pi);

            ug =  u_loc*cb - v_loc*sb;
            vg =  u_loc*sb + v_loc*cb;

            Uind_t = Uind_t + (ug*tx(i) + vg*ty(i))*sigma(j);
        end
        Ut(i) = U_t + Uind_t;
    end

    Cp = 1 - (Ut./Uinf).^2;
end