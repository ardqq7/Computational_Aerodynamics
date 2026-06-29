function [Vt, Cp] = panelVelocityCp(Uinf, alphaDeg, geom, Kij, lambda)
% Compute tangential surface velocity and Cp using Kij and lambda_j

    N = geom.N;

    % Free-stream direction: here we assume Vinf along x
    alphaRad = deg2rad(alphaDeg);
    Ufx = Uinf * cos(alphaRad);
    Ufy = Uinf * sin(alphaRad);

    % Tangent unit vector at each collocation point
    tx = geom.tx(:);
    ty = geom.ty(:);

    % Free-stream tangential component
    Vt_inf = Ufx .* tx + Ufy .* ty;   % N x 1

    % Panel strengths (if lambda is N x N influence, you probably have a
    % solved vector lambda_j somewhere; here I assume you already solved
    % the system and have lambda_j as N x 1)
    %
    % If your notation is: Vij = (1/(2π)) * sum_j Kij * λ_j + V_inf,
    % then:
    lambda_j = diag(lambda);  % <-- örnek: eğer lambda_ij değil de λ_j kullanıyorsan
                              % kendi sistemine göre değiştir.

    % Influence of panels on tangential velocity
    Vt_induced = (1/(2*pi)) * (Kij * lambda_j);   % N x 1

    % Total tangential velocity
    Vt = Vt_inf + Vt_induced;

    % Pressure coefficient from Bernoulli:
    Cp = 1 - (Vt./Uinf).^2;
end