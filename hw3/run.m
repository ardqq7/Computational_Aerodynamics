clear; clc; close all;

% Geometry parameters
r   = 10;
cx  = 0;
cy  = 0;
startAngleDeg = 157.5;

% Flow conditions
Uinf     = 1.0;
alphaDeg = 0.0;   % saldırı açısı

% Çözmek istediğin panel sayıları
N_list = [8 32 64];

% Sonuçları saklamak için
Cp_all     = cell(numel(N_list),1);
theta_all  = cell(numel(N_list),1);

for k = 1:numel(N_list)
    N = N_list(k);

    % --- Panelleri üret (senin fonksiyonun) ---
    fileName = sprintf('circle_points_N%d.csv', N);
    [x, y, ~, ~, ~, ~, ~, ~, ~] = makeCirclePanels(r, N, cx, cy, startAngleDeg, fileName);

    % Emin ol x,y kapalı contour: son nokta = ilk nokta
    % if x(1) ~= x(end) || y(1) ~= y(end)
    %     x = [x(:); x(1)];
    %     y = [y(:); y(1)];
    % else
    %     x = x(:); y = y(:);
    % end

    % --- Geometriyi kur ---
    % buildPanelGeometry'yi şöyle tanımla:
    % [geom, xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);
    [geom, xc, yc, S, beta, phi, tx, ty, nx, ny] = buildPanelGeometry(x, y);
    % geom zaten bunları içeren struct olsun:
    % geom.x, geom.y, geom.midpp, geom.S, geom.tx, geom.ty, geom.nx, geom.ny, geom.beta, geom.phi, geom.N

    % --- Etkileşim katsayıları (Katz & Plotkin 3.158–3.163) ---
    ABCDES = buildABCDES(geom);
    [Iij, Aij] = buildI_and_A(ABCDES); %#ok<NASGU> % Iij istersen sonra kullanırsın

    % --- RHS: no-penetration: U · n + sum(sigma_j * I_ij) = 0 ---
    Ufx = Uinf*cosd(alphaDeg);
    Ufy = Uinf*sind(alphaDeg);
    b   = -(Ufx*geom.nx + Ufy*geom.ny);   % N×1

    % Kapanış: sum(sigma)=0
    Aij(end,:) = 1.0;
    b(end)     = 0.0;

    % --- Kaynak şiddetleri ---
    sigma = Aij \ b;

    % --- Teğetsel hız ve Cp ---
    [Ut, Cp] = panelTangentialCp(Uinf, alphaDeg, geom, sigma);

    % --- Her panel orta noktası için kutupsal açı θ ---
    xc_mid = geom.midpp(:,1);
    yc_mid = geom.midpp(:,2);
    theta = atan2(yc_mid - cy, xc_mid - cx);   % [-pi,pi]

    % Sonuçları sakla
    Cp_all{k}    = Cp;
    theta_all{k} = theta;
end

%% 2) Teorik Cp (kaldırmasız silindir, potansiyel akım)
theta_th = linspace(0, 2*pi, 400);
alphaRad = deg2rad(alphaDeg);
Cp_th = 1 - 4*(sin(theta_th - alphaRad)).^2;

%% 3) Hepsini aynı grafikte çiz
figure(10); clf;
set(gcf,'Color','w'); hold on; grid on; box on;

% Teorik
plot(rad2deg(theta_th), Cp_th, 'k-', 'LineWidth', 1.8, 'DisplayName','Theory');

% Sayısal (her N için)
clr = lines(numel(N_list));
for k = 1:numel(N_list)
    th_deg = rad2deg(theta_all{k});
    
    % Açıyı [0,360) aralığına al
    th_deg = mod(th_deg,360);
    
    % Sıralayalım ki grafik düzgün olsun
    [th_sorted, idx] = sort(th_deg);
    Cp_sorted = Cp_all{k}(idx);
    
    plot(th_sorted, Cp_sorted, 'o-','LineWidth',1.3, ...
        'Color',clr(k,:), ...
        'MarkerSize',4, ...
        'DisplayName', sprintf('N = %d', N_list(k)));
end

xlabel('\theta (deg)');
ylabel('C_p');
title('C_p around a circular cylinder (source panel method)');
set(gca,'XLim',[0 360]);
legend('Location','best');