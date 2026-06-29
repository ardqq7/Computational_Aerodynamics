%% Arda Ceker 110200131
% UCK 419E - Computational Aerodynamics - Term Project - 2025/2026 Fall
clear, clc, clf, close all, format long g;

%% Physical parameters from ISA at sea level
[T,a,P,rho,nu,mu] = atmosisa(0); % ISA, sea level
U_inf = 1.0; % Free stream velocity in m/s

%% Domain size
Lx = 1.0; % domain length in x-direction in m
Ly = 0.05; % domain height in y-direction in m

%% Grid (pressure control volumes)
Nx = 200; % number of cells in x
Ny = 200; % number of cells in y
    
dx = Lx / Nx;
dy = Ly / Ny;

Re_L = U_inf * Lx / nu; % Reynolds number

%% Time parameters
dt = 1.0e-4; % time step in s
t_end = 1.5; % final time ins
nt = round(t_end / dt);

%% Boundary conditions
U_inlet = U_inf; % inlet 
V_inlet = 0.0; % inlet 
U_top = U_inf; % free-stream
U_bottom = 0.0; % no-slip wall

%% Coordinates
% Pressure at cell centers
xp = linspace(dx/2, Lx - dx/2, Nx); % Cell centers
yp = linspace(dy/2, Ly - dy/2, Ny); % Cell centers

% u-velocity at vertical faces
xu = linspace(0, Lx, Nx+1);
yu = yp; % same y as pressure centers

% v-velocity at horizontal faces
xv = xp; % same x as pressure centers
yv = linspace(0, Ly, Ny+1);

p = zeros(Nx,   Ny); % pressure at cell centers
u = zeros(Nx+1, Ny); % x-velocity at vertical faces
v = zeros(Nx,   Ny+1); % y-velocity at horizontal faces

%% Initial conditions
u(:) = U_inf;
v(:) = 0.0;
p(:) = 0.0;

u(1, :) = U_inlet; % Free stream BC
v(1, :) = V_inlet;  % Free stream BC
v(:, 1) = 0.0; % No slip
u(:, 1) = U_bottom; % No slip
v(:, Ny+1)= 0.0; % No slip
u(:, Ny) = U_top; % Free stream BC

u_star = u;
v_star = v;
p_corr = zeros(Nx, Ny);

u_old = u;
v_old = v;
p_old = p;

%% Analysis Parameters
fprintf('Re_L = %.2f\n', Re_L);
fprintf('dx, dy = %.4e, %.4e\n', dx, dy);
fprintf('dt = %.1e, t_end = %.2f, nt = %d\n', dt, t_end, nt);

%% Solution
tic; % Computation time is start

for n = 1:nt
    u(1, :)  = U_inlet; % inlet
    u(:, 1)  = U_bottom; % wall
    u(:, Ny) = U_top; % top (free-stream)
    u_old = u;
    v(:,1) = 0.0;

    for i = 1:Nx % velocity u-cells
        for j = 1:Ny % v-cells
            du_dx = (u_old(i+1,j) - u_old(i,j)) / dx; % velocity change
            v(i,j+1) = v(i,j) - dy * du_dx; % next cell calculation for v
        end
    end

    % u update
    u_new = u_old;

    for i = 2:Nx % interior u-faces
        for j = 2:Ny-1 % v-faces
            u_ij = u_old(i,j);
            du_dx = (u_old(i+1,j) - u_old(i-1,j)) / (2*dx); % velocity change
            du_dy = (u_old(i,j+1) - u_old(i,j-1)) / (2*dy); % velocity change
            v_u = 0.5 * ( v(i-1,j) + v(i,j) );  % v at u location 

            u_yy = (u_old(i,j+1) - 2*u_old(i,j) + u_old(i,j-1)) / dy^2; 

            conv = u_ij * du_dx + v_u * du_dy;
            diff = nu * u_yy;

            u_new(i,j) = u_ij + dt * ( -conv + diff );

        end
    end

    % Boundary conditions
    u_new(1, :)  = U_inlet;
    u_new(:, 1)  = U_bottom;
    u_new(:, Ny) = U_top;
    u = u_new;

    if mod(n, floor(nt/10)) == 0
        mid_i = round(Nx/2);
        figure(1); clf; % Time dependent velocity profiles
        plot(u(mid_i,:), yu, 'LineWidth', 2);
        grid minor;
        xlabel('u [m/s]');
        ylabel('y [m]');
        title(sprintf('Velocity Profile t = %.4f s, at x = %.3f m', n*dt, xu(mid_i)));
        ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
        set(gcf, 'Color' , 'w')
        xlim([0 1.2])
        drawnow;
    end

end

elapsed_time = toc;  % Computing time in total s
fprintf('Compute time = %.4f s\n', elapsed_time);

%% Results

u_c = 0.5 * (u(1:Nx,:) + u(2:Nx+1,:)); % Veloctiy u at the center cell

delta_num = zeros(1, Nx); % Data store
delta_star_num = zeros(1, Nx); % Data store
theta_num = zeros(1, Nx); % Data store
Cf_num = zeros(1, Nx); % Data store

for i = 1:Nx
    u_col = u_c(i,:) / U_inf; % Velocity porofile at center cell
    idx = find(u_col >= 0.99, 1, 'first'); % Boundary layer at 0.99u_inf
    if isempty(idx)
        delta_num(i) = Ly;
    else
        delta_num(i) = yu(idx);
    end

    delta_star_num(i) = sum( (1 - u_col) ) * dy;  % Integrateing displacement thickness
    theta_num(i) = sum( u_col .* (1 - u_col) ) * dy;  % Integrateing momentunm thickness
    % du/dy at wall
    du_dy_wall = (u_c(i,2) - u_c(i,1)) / dy; % Shear stress
    tau_w      = mu * du_dy_wall; % % Shear stress
    Cf_num(i)  = 2 * tau_w / (rho * U_inf^2); % Skin friction coefficient
end

%% Boundary Layer Parameters by Theoretical / Analytical Approach - Blasius Solution
Re_x = U_inf .* xp / nu; % Reynolds number
delta_th = 5   * sqrt(nu * xp / U_inf); % Boundary layer thickness
deltaS_th = 1.72* sqrt(nu * xp / U_inf); % Displacement thickness
theta_th = 0.664*sqrt(nu * xp / U_inf); % Momentum thickness
Cf_th = 0.664 ./ sqrt(Re_x); % Skin friction coefficient

%% Plot: delta(x)
figure(2);
gca;
plot(xp, delta_num, 'LineWidth', 2); hold on;
plot(xp, delta_th,  '--', 'LineWidth', 2);
xlabel('x [m]');
ylabel('Boundary Layer Thickness \delta [m]');
legend('Numerical Solution','Analytical Solution','Location','best');
title('Boundary Layer Thickness \delta(x) along the Plate');
grid minor;
ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
lgd = legend('Location','best'); lgd.FontWeight = 'bold'; lgd.FontSize = 11;
set(gcf, 'Color' , 'w')

%% Displacement thickness delta-star - momentum thickness theta

figure(3);
    plot(xp, delta_star_num, 'LineWidth', 2); hold on;
    plot(xp, theta_num,      'LineWidth', 2);
    plot(xp, deltaS_th, '--', 'LineWidth', 1.5);
    plot(xp, theta_th,  '--', 'LineWidth', 1.5);
    grid minor;
    xlabel('x [m]');
    ylabel('Thickness [m]');
    legend('\delta^*_{numerical}','\theta_{numerical}', ...
           '\delta^*_{analytical}','\theta_{analytical}','Location','best');
    title('Displacement Thickness \delta^*(x) and Momentum Thickness \theta(x)');
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    lgd = legend('Location','best'); lgd.FontWeight = 'bold'; lgd.FontSize = 11;
    set(gcf, 'Color' , 'w')

%% Skin friction coefficient Cf

figure(4);
    plot(xp, Cf_num, 'LineWidth', 2); hold on;
    plot(xp, Cf_th,  '--', 'LineWidth', 2);
    grid on;
    xlabel('x [m]');
    ylabel('C_f');
    legend('Numerical Solution','Analytical Solution','Location','best');
    title('Skin Friction Coefficient C_f(x)');
    grid minor;
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    lgd = legend('Location','best'); lgd.FontWeight = 'bold'; lgd.FontSize = 11;
    set(gcf, 'Color' , 'w')
    grid minor;
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    lgd = legend('Location','best'); lgd.FontWeight = 'bold'; lgd.FontSize = 11;
    set(gcf, 'Color' , 'w')

%% Contour of u/U_inf
figure(5);
    [XX,YY] = meshgrid(xp, yu);
    contourf(XX,YY,(u_c'/U_inf),20,'LineStyle','none');
    colormap turbo
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    title('u/U_\infty contour');
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    set(gcf, 'Color' , 'w')

%% Pressure contour 

[XXc, YYc] = meshgrid(xp, yp); % pressure at cell centersi in PA

figure(6);
    contourf(XXc, YYc, p', 20, 'LineStyle', 'none');
    colormap turbo
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    title('Pressure contour');
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    set(gcf, 'Color' , 'w')

%% Non-dimensional velocity profiles

i_profile = round(Nx/2); % Center of the plate
x_loc = xp(i_profile); % Position of the center in m
u_c = 0.5 * (u(1:Nx,:) + u(2:Nx+1,:)); % Velocity profiu at the center
u_num = u_c(i_profile,:) / U_inf; % Non dimensional velocity u/U_inf

delta_loc_num = delta_num(i_profile); % Boundary layer thickness at 99% u_inf grid no
y_num = yu; % non-dimensional distance from wall
eta_num = y_num / delta_loc_num; % eta = y/delta

mask = eta_num <= 1.0; % Mask for 0 <= y/delta
eta_num_plot = eta_num(mask); % Only masked elements
u_num_plot = u_num(mask);

eta_th = linspace(0,1,200); % Non dimensional eta
u_lam = 1.5*eta_th - 0.5*eta_th.^3; % Laminar approximation
u_turb = eta_th.^(1/7); % Turbulent approximation

figure(7); % Plotting
    plot(u_num_plot, eta_num_plot, 'o-','LineWidth',1.5); hold on;
    plot(u_lam,  eta_th, '-','LineWidth',2);
    plot(u_turb, eta_th, '--','LineWidth',2);
    grid minor;
    axis equal;
    xlabel('u/U_\infty');
    ylabel('y/\delta');
    title(sprintf('Non-dimensional Velocity Profile at x = %.3f m', x_loc));
    xlim([0 1]); ylim([0 1]); 
    ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 12; ax.YAxis.Exponent = 0;
    set(gcf, 'Color' , 'w')
    xlim([0 1]); ylim([0 1]);
    legend('Numerical Solution', 'Laminar Boundary Layer', 'Turbulent Boundary Layer', 'Location','best');

%% Integral mass and momentum transfer through the boundaries

m_dot_inlet = rho * sum( u(1,:) ) * dy;
m_dot_exit = rho * sum( u(Nx+1,:) ) * dy;
m_dot_bottom = rho * sum( v(:,1) ) * dx;
m_dot_top = rho * sum( v(:,Ny+1) ) * dx;
m_dot_net = m_dot_inlet + m_dot_bottom - m_dot_exit - m_dot_top; % Mass flow rate in m/s

Mx_inlet  = rho * sum( (u(1,:).^2)    ) * dy;
Mx_exit   = rho * sum( (u(Nx+1,:).^2) ) * dy;
Mx_bottom = rho * sum( u_c(:,1) .* v(:,1) ) * dx;
Mx_top    = rho * sum( u_c(:,Ny) .* v(:,Ny+1) ) * dx;

%% Momentum transfer to the wall shear
u_wall_grad = (u_c(:,2) - u_c(:,1)) / dy;
tau_w_dist  = mu * u_wall_grad; % wall shear stress 
D_wall = sum(tau_w_dist) * dx; % drag force on the wall

%% Print integral quantities
fprintf('\nStreamwise mass flux :\n');
fprintf('Mass flux at inlet : % .4e kg/s per unit span\n', m_dot_inlet);
fprintf('Mass flux at outlet : % .4e kg/s per unit span\n', m_dot_exit);
fprintf('Net mass flow between boundaries : % .4e kg/s per unit span\n', m_dot_net);

fprintf('\nStreamwise momentum flux :\n');
fprintf('x-momentum at inlet : % .4e N per unit span\n', Mx_inlet);
fprintf('x-momentum at outlet : % .4e N per unit span\n', Mx_exit);
fprintf('x-momentum at bottom : % .4e N per unit span\n', Mx_bottom);
fprintf('x-momentum at top : % .4e N per unit span\n', Mx_top);

fprintf('\nMomentum transfer to the wall:\n');
fprintf('Total wall shear force (drag) : % .4e N per unit span\n', D_wall);

%% Solution Stability
CFL_x_Uinf = U_inf * dt / dx; % convection
CFL_y_Uinf = U_inf * dt / dy; % diffusion
Fo_y       = nu    * dt / dy^2; % diffusion only in y-direction

u_max = max(abs(u(:))); % maximum CFLx
v_max = max(abs(v(:))); % maximum CFLy
CFL_x_max = u_max * dt / dx;
CFL_y_max = v_max * dt / dy;

fprintf('CFL_x = %.3e\n', CFL_x_Uinf);
fprintf('CFL_y = %.3e\n', CFL_y_Uinf);
fprintf('Fo_y = %.3e\n', Fo_y);

fprintf('CFL_x(max|u|) = %.3e\n', CFL_x_max);
fprintf('CFL_y(max|v|) = %.3e\n', CFL_y_max);

%% Error of Boundary Layer Thickness - Analytical and Numerical
i_mid = round(Nx/2);
x_mid = xp(i_mid); % Domain center

Re_x_mid = U_inf * x_mid / nu;
delta_mid_num = delta_num(i_mid); % numerical solution (delta_99)
delta_mid_th = delta_th(i_mid);  % Blasius 

rel_err_delta = (delta_mid_num - delta_mid_th) / delta_mid_th * 100; % error

fprintf('\nBoundary-layer thickness at mid-domain (x = %.3f m)\n', x_mid);
fprintf('Numerical = %.4e m\n', delta_mid_num);
fprintf('Blasius = %.4e m\n', delta_mid_th);
fprintf('Relative error (numerical vs. Blasius) = %.2f %%\n', rel_err_delta);
