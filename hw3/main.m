% Arda Çeker 110200131
% Homework 3 - Code 

clear, clc, clf, close all, format long g;

% Panel sayıları
N_list = [8 64 256];

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

    % panelcode fonksiyonunun N için Cp ve theta üretmesi gerekir
    [Cp, theta_deg] = panelcode(N);  

    % Sort for smooth plot
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

%% --- Analitik çözüm ---
theta_th = linspace(0,360,400);
theta_rad = deg2rad(theta_th);

% alpha = 0 için bilinen formül
Cp_th = 1 - 4*(sin(theta_rad)).^2;
figure(100)
plot(theta_th, Cp_th, 'k-', 'LineWidth', 1.8, 'DisplayName','Analytical');

legend('Location','best');
xlim([0 360]);