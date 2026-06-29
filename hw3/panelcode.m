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

    % Close contour if needed
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

    % Optional: label panel indices at midpoints
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

% Geometriyi çıkar
[geom, xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);

% midpp erişimi:
midpp = geom.midpp;     % Nx2, [xc yc]
% İstersen ayrı ayrı:
xc = midpp(:,1); 
yc = midpp(:,2);

% Çizdir (teğet/normal oklarıyla)
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

% (D) Hızlı grafik (opsiyonel)
figure(4);
set(gcf, 'Color','w'); plot(1:geom.N, Cpp, '-o','LineWidth',1.5); grid on;
xlabel('Panel index'); ylabel('C_p'); title('Surface C_p (source panels)');

% Geometri (x,y kapalı; CCW önerilir)
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

% Kapanış şartı: sum(sigma)=0 → son satırı bununla değiştir
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

