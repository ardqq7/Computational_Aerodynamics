function ANG = buildAngleMatrix(geom)
% Build angle-related matrices (beta_i, phi_j, beta_i - phi_j, etc.)

    N = geom.N;

    beta = geom.beta(:);  % tangent angle at control point i
    phi  = geom.phi(:);   % normal angle of panel j

    % Repeat to form NxN matrices
    BETA_i = repmat(beta, 1, N);   % each column = beta_i
    PHI_j  = repmat(phi.', N, 1);  % each row    = phi_j

    % Angle differences you'll probably need
    dBetaPhi = BETA_i - PHI_j;     % β_i - Φ_j

    ANG = struct();
    ANG.BETA_i  = BETA_i;
    ANG.PHI_j   = PHI_j;
    ANG.dBetaPhi = dBetaPhi;
end