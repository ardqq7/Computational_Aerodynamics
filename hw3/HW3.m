% Arda Çeker 110200131
% Homework 3 - Code 

clear, clc, clf, close all, format long g;

N_list = [8 32 64 256]; % Panels

% Circle center
cx = 0; 
cy = 0;

figure; hold on;
set(gcf,'Color','w');
grid on; box on;
xlabel('\theta (deg)');
ylabel('C_p');
title('Cp(\theta) for different panel counts');

colors = lines(numel(N_list));

for k = 1:numel(N_list)

    N = N_list(k);

    [Cp, theta_deg] = panelcode(N);  

    % [theta_sorted, idx] = sort(theta_deg);
    % Cp_sorted = Cp(idx);
    figure(100);
    hold on
    plot(theta_deg, Cp, '-o', ...
         'LineWidth', 1.2, 'MarkerSize', 4, ...
         'Color', colors(k,:), ...
         'DisplayName', sprintf('N = %d', N));


    set(gcf,'Color','w');
    % plot(theta_deg, Cp, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    grid on; box on;     
    xlabel('\theta (deg)');
    ylabel('C_p');
    title('C_p distribution around the circle');
    xlim([0 360]);
end

% --- Analytical ---
theta_th = linspace(0,360,400);
theta_rad = deg2rad(theta_th);

Cp_th = 1 - 4*(sin(theta_rad)).^2;
figure(100)
plot(theta_th, Cp_th, 'k-', 'LineWidth', 1.8, 'DisplayName','Analytical');

legend('Location','best');
xlim([0 360]);

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

function plotPanelsMidTangNorm(geom, scaleTN)

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

function [sigma, A, b] = assembleSourceSystem(Uinf, alphaDeg, geom)
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

function ABCDES = buildABCDES(geom)
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

function [Iij, Aij] = buildI_and_A(ABCDES)

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

function [Ut, Cp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma)
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

function [Cpp, theta_deg] = panelcode(N)
% clear, clc, clf, close all, format long g;

%data = readmatrix('input.txt'); % Reads file and saves data as a matrix
%x = data(:,1); % Stores x values from first column
%y = data(:,2); % Stores y values from first column

r = 10;

cx = 0;
cy = 0;
startAngleDeg = 180;
fileName = 'circle_points.csv';

[x, y, xc, yc, tx, ty, nx, ny, S] = makeCirclePanels(r, N, cx, cy, startAngleDeg, fileName);

x = x(:); y = y(:);

    % Close contour
    %if x(1)~=x(end) || y(1)~=y(end)
    %    x = [x; x(1)];
    %    y = [y; y(1)];
   % end

    % Panel count
    N = numel(x)-1;

    % Midpoints (collocation points)
    xc = (x(1:end-1) + x(2:end)) / 2;
    yc = (y(1:end-1) + y(2:end)) / 2;

    % --- Plot ---
    figure(1);
    set(gcf, 'Color','w'); hold on; axis equal; grid on; box on;

    % Panels (connect nodes in order)
    plot(x, y, 'k-', 'LineWidth', 1.5);

    % Nodes
    plot(x(1:end-1), y(1:end-1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

    % Midpoints in red
    plot(xc, yc, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

    for i = 1:N
        text(xc(i), yc(i), sprintf('  %d', i), 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    end

    % Reference circle (green dashed)
    th = linspace(0, 2*pi, 400);
    xcirc = cx + r*cos(th);
    ycirc = cy + r*sin(th);
    plot(xcirc, ycirc, 'g--', 'LineWidth', 1.2);

   
    title(sprintf('Panels and Midpoints (N = %d)', N));
    xlabel('x'); ylabel('y');
    legend({'Panels','Nodes','Midpoints','Reference Circle'}, 'Location','best');

[x,y,~,~,~,~,~,~,~] = makeCirclePanels(r, N, cx, cy, startAngleDeg, fileName);  % senin fonksiyonun

% Geometry
[geom, xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);

% midp
midpp = geom.midpp;     % Nx2, [xc yc]

xc = midpp(:,1); 
yc = midpp(:,2);

plotPanelsMidTangNorm(geom, 0.35);

    beta_deg = rad2deg(beta(:));
    phi_deg  = rad2deg(phi(:));

    % Unwrap to avoid jumps at +/-180 deg
    beta_deg_u = unwrap(deg2rad(beta_deg))*180/pi;
    phi_deg_u  = unwrap(deg2rad(phi_deg))*180/pi;

    N = numel(beta_deg);

    figure(3); clf; set(gcf,'Color','w'); hold on; grid on; box on;
    plot(1:N, beta_deg_u, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    plot(1:N,  phi_deg_u, '-s', 'LineWidth', 1.5, 'MarkerSize', 4);

    xlabel('Panel index');
    ylabel('Angle (deg)');
    title('\beta (tangent) and \phi (normal) angles vs panel index');
    legend({'\beta (deg)','\phi (deg)'}, 'Location', 'best');

    % (A) Circle or arbitrary shape already discretized:
% x, y from your generator or file (closed, CCW recommended)
[xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);
%geom = struct('x',x,'y',y,'midpp',[xc yc],'S',S,'tx',tx,'ty',ty,'nx',nx,'ny',ny,'beta',beta,'phi',phi,'N',numel(S));

% (B) Assemble and solve for sigma
Uinf = 1.0; alphaDeg = 0.0;
[sigma, A, b] = assembleSourceSystem(Uinf, alphaDeg, geom);

% (C) Tangential velocity and Cp at midpoints
[Ut, Cpp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma);

% (D) 
figure(4);
set(gcf, 'Color','w'); plot(1:geom.N, Cpp, '-o','LineWidth',1.5); grid on;
xlabel('Panel index'); ylabel('C_p'); title('Surface C_p (source panels)');

% Geometry
[xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);
%geom = struct('x',x,'y',y,'midpp',[xc yc],'S',S,'tx',tx,'ty',ty,...
      %        'nx',nx,'ny',ny,'beta',beta,'phi',phi,'N',numel(S));

ABCDES = buildABCDES(geom);

% I_ij ve A_ij (A_ii=1/2 olacak şekilde)
[Iij, Aij] = buildI_and_A(ABCDES);

% RHS b_i = -U∞ · n_i
Uinf = 10.0; alphaDeg = 0.0;
Ufx = Uinf*cosd(alphaDeg); Ufy = Uinf*sind(alphaDeg);
b = -(Ufx*nx + Ufy*ny);

% Kapanış şartı: sum(sigma)=0
Aij(end,:) = 1.0; b(end) = 0.0;

sigma = Aij \ b;   % Source

        % [Ut, Cp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma);
ANG = buildAngleMatrix(geom);

[lambdaMat, Kij] = buildLambdaAndK(Aij, ANG);

lambda_j = sigma;  

[Vt, Cp] = panelVelocityCp(Uinf, alphaDeg, geom, Kij, diag(lambda_j));

% Center of the circle
cx = 0;   
cy = 0;
[Ut, Cp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma);

% Midpoints of panels
xc_mid = geom.midpp(:,1);
yc_mid = geom.midpp(:,2);

% Polar angle of each midpoint (radians)
theta = atan2(yc_mid - cy, xc_mid - cx);   % [-pi, pi]
theta = theta + pi;
% Convert to degrees and shift to [0, 360)
theta_deg = rad2deg(theta);
% theta_deg = mod(theta_deg, 360);

% Sort by theta so that the curve is smooth
[theta_sorted, idx] = sort(theta_deg);
Cp_sorted = Cp(idx);

test = 0;
% Plot Cp vs theta
% figure(5);
% set(gcf, 'Color','w');
% plot(theta_sorted, Cp_sorted, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
% grid on; box on;
% xlabel('\theta (deg)');
% ylabel('C_p');
% title('C_p distribution around the circle');
% xlim([0 360]);
end

function ANG = buildAngleMatrix(geom)
% Build angle-related matrices (beta_i, phi_j, beta_i - phi_j, etc.)

    N = geom.N;

    beta = geom.beta(:);  % tangent angle at control point i
    phi  = geom.phi(:);   % normal angle of panel j

    % Repeat to form NxN matrices
    BETA_i = repmat(beta, 1, N);   % each column = beta_i
    PHI_j  = repmat(phi.', N, 1);  % each row    = phi_j

    % Angle differences 
    dBetaPhi = BETA_i - PHI_j;     % β_i - Φ_j

    ANG = struct();
    ANG.BETA_i  = BETA_i;
    ANG.PHI_j   = PHI_j;
    ANG.dBetaPhi = dBetaPhi;
end

function [lambda, Kij] = buildLambdaAndK(Aij, ANG)
% From normal influence Aij and angle matrix, build lambda_ij and Kij

    dBetaPhi = ANG.dBetaPhi;
    N = size(Aij,1);

    lambda = zeros(N);   % lambda_ij
    Kij    = zeros(N);   % K_ij

    % lambda_ij = A_ij * cos(β_i - Φ_j)
    lambda = Aij .* cos(dBetaPhi);

    % K_ij = A_ij * sin(β_i - Φ_j) 
    Kij = Aij .* sin(dBetaPhi);
    % ---------------------------------------------------

end

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

    % Panel strengths
    %  Vij = (1/(2π)) * sum_j Kij * λ_j + V_inf,

    lambda_j = diag(lambda); 

    % Influence of panels on tangential velocity
    Vt_induced = (1/(2*pi)) * (Kij * lambda_j);   % N x 1

    % Total tangential velocity
    Vt = Vt_inf + Vt_induced;

    % Pressure coefficient from Bernoulli:
    Cp = 1 - (Vt./Uinf).^2;
end