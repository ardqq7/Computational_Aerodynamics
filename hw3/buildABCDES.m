function ABCDES = buildABCDES(geom)
% Build A,B,C,D,E,S matrices for Eq. (3.162)–(3.163)
% geom fields used: midpp(:,1/2)=xi,yi ; x,y ; phi (Nx1); S (Nx1)
% Notation:
%   (xi,yi)  : i-th control (collocation) point
%   (Xj,Yj)  : j-th panel START node
%   Sj       : j-th panel length
%   phi_i    : normal angle at i (rad)
%   phi_j    : normal angle of panel j (rad)

    xi  = geom.midpp(:,1);   yi  = geom.midpp(:,2);   % size N
    Xj  = geom.x(1:end-1);   Yj  = geom.y(1:end-1);   % size N (panel starts)
    Sj  = geom.S(:);                               % size N
    phi = geom.phi(:);                             % size N (normal angle)

    N = geom.N;
    % Preallocate
    A = zeros(N); B = zeros(N); C = zeros(N); D = zeros(N);
    E = zeros(N); S = repmat(Sj.', N, 1);          % broadcast Sj along rows

    % Build pairwise differences Δx, Δy between (xi,yi) and panel-j start (Xj,Yj)
    DX = xi - Xj.';     % N x N, i rows, j cols
    DY = yi - Yj.';     % N x N

    % Make per-i and per-j angle matrices
    PHI_i = repmat(phi, 1, N);      % each column copies phi_i
    PHI_j = repmat(phi.', N, 1);    % each row    copies phi_j

    % ----- Formulas from the screenshot -----
    % B = (xi - Xj)^2 + (yi - Yj)^2
    B = DX.^2 + DY.^2;

    % A = -(xi - Xj) cos(phi_j) - (yi - Yj) sin(phi_j)
    A = -(DX .* cos(PHI_j) + DY .* sin(PHI_j));

    % C = sin(phi_i - phi_j)
    C = sin(PHI_i - PHI_j);

    % D = (yi - Yj) cos(phi_i) - (xi - Xj) sin(phi_i)
    D = DY .* cos(PHI_i) - DX .* sin(PHI_i);

    % E = sqrt(B - A^2)  (≡ (xi - Xj) sin φ_j − (yi − Yj) cos φ_j up to sign)
    E = sqrt(max(B - A.^2, 0));     % numerical safety

    ABCDES = struct('A',A,'B',B,'C',C,'D',D,'E',E,'S',S,'DX',DX,'DY',DY);
end