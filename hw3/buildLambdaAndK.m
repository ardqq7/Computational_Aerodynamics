function [lambda, Kij] = buildLambdaAndK(Aij, ANG)
% From normal influence Aij and angle matrix, build lambda_ij and Kij
%
% NOTE:
%   Replace the placeholder formulas with the exact ones from your notes.

    dBetaPhi = ANG.dBetaPhi;
    N = size(Aij,1);

    lambda = zeros(N);   % lambda_ij
    Kij    = zeros(N);   % K_ij

    % ----- PLACEHOLDER: put your own formulas here -----
    % ÖRNEK (tamamen uydurma, SADECE iskelet):
    % lambda_ij = A_ij * cos(β_i - Φ_j)
    lambda = Aij .* cos(dBetaPhi);

    % K_ij için de benzer şekilde:
    % K_ij = A_ij * sin(β_i - Φ_j)   (örnek!)
    Kij = Aij .* sin(dBetaPhi);
    % ---------------------------------------------------

end