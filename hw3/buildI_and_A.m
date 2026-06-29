function [Iij, Aij] = buildI_and_A(ABCDES)
% Compute I_ij via Eq. (3.163) and A = I/(2π) with A_ii = 1/2

    A = ABCDES.A; B = ABCDES.B; C = ABCDES.C; D = ABCDES.D;
    E = ABCDES.E; S = ABCDES.S;          % S is N x N with column j = Sj

    N = size(A,1);
    Iij = zeros(N);

    % Avoid division-zero; use limiting forms when E ~ 0
    epsE = 1e-14; E_safe = E; E_safe(E_safe < epsE) = epsE;

    % Eq. (3.163):
    % I_ij = (C/2) * ln( (S_j^2 + 2 A S_j + B)/B )  + (D - A C)/E * ( atan((S_j + A)/E) - atan(A/E) )
    termLog = log( (S.^2 + 2.*A.*S + B) ./ B );
    termAtn = atan( (S + A) ./ E_safe ) - atan( A ./ E_safe );

    Iij = 0.5.*C .* termLog + ( (D - A.*C) ./ E_safe ) .* termAtn;

    % Influence matrix for no-penetration using source panels:
    Aij = Iij ./ (2*pi);

    % Diagonal (self-influence) limit -> 1/2
    for i = 1:N
        Aij(i,i) = 0.5;
    end
end