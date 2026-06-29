%% Arda Çeker 110200131
% UCK419E - Computational Aerodynamics - HW4 

clear; clc; close all; format long g;

%% 1) Physical properties - Aluminum
rho   = 2700;        % [kg/m^3] density
cp    = 900;         % [J/kg-K] specific heat
k     = 230;        % [W/m-K] thermal conductivity
alpha = k/(rho*cp);  % [m^2/s] thermal diffusivity (just for info)

%% 2) Geometry and BC parameters (common for all runs)
Lx = 0.5;            % [m] plate length in x-direction
Ly = 0.5;            % [m] plate length in y-direction

T_initial = 300.0;   % [K] initial temperature in the plate
T_bottom  = 400.0;   % [K] bottom edge (y = 0) - heated surface
T_top     = 300.0;   % [K] top edge (y = Ly)
T_left    = 300.0;   % [K] left edge (x = 0)
T_right   = 300.0;   % [K] right edge (x = Lx)

% param struct: input for the function
params.rho       = rho;
params.cp        = cp;
params.k         = k;
params.Lx        = Lx;
params.Ly        = Ly;
params.T_initial = T_initial;
params.T_bottom  = T_bottom;
params.T_top     = T_top;
params.T_left    = T_left;
params.T_right   = T_right;

%% 3) Base case for showing results
Nx0      = 41;       % base number of nodes in x
Ny0      = 41;       % base number of nodes in y
dt0      = 0.1;     % base time step [s]
t_final  = 100.0;     % total simulation time [s]

[T_final_base, x_base, y_base, t_hist_base, T_center_hist_base, Q_bottom_hist_base, cpu_time_base] = HeatTransferAnalysis(Nx0, Ny0, dt0, t_final, params, true); % calling anaylsis function

fprintf('\nResults for : Nx = %d, Ny = %d, dt = %.3f, CPU time = %.3f s\n', Nx0, Ny0, dt0, cpu_time_base);


%% 4) Grid & time-step study for Mesh and Time independency 
Nx_ref = 101; % reference grid in x
Ny_ref = 101; % reference grid in y
dt_ref = 0.005; % reference time step [s]

[~, ~, ~, t_hist_ref, T_center_hist_ref, ~, cpu_time_ref] = HeatTransferAnalysis(Nx_ref, Ny_ref, dt_ref, t_final, params, true); % calling anaylsis function

T_center_ref_final = T_center_hist_ref(end); % Final temperature
fprintf('Base Computation Time: Nx = %d, Ny = %d, dt = %.3f, CPU time = %.4f s\n', Nx_ref, Ny_ref, dt_ref, cpu_time_ref);

Nx_list = [21, 31, 41, 51];      % Nx = Ny square grid
dt_list = [0.10, 0.05, 0.025];   % Time step comparison

% Results matrices
nNx = numel(Nx_list);
nDt = numel(dt_list);

err_center = zeros(nNx, nDt);   % Center temperature error at t_final
cpu_time   = zeros(nNx, nDt);   % CPU time for each run

for ii = 1:nNx
    for jj = 1:nDt

        Nx_test = Nx_list(ii);
        Ny_test = Nx_list(ii);   % Nx = Ny
        dt_test = dt_list(jj);

        [~, ~, ~, ~, T_center_hist_test, ~, cpu_time_test] = HeatTransferAnalysis(Nx_test, Ny_test, dt_test, t_final, params, false); % calling anaylsis function

        T_center_test_final = T_center_hist_test(end); % Final temperature

        err_center(ii,jj) = abs(T_center_test_final - T_center_ref_final); % Error between reference (fine) and this one
        cpu_time(ii,jj)   = cpu_time_test; % Elapsted time

        fprintf(' -> T_center(t_final) = %.6f K, error = %.6f K, CPU = %.4f s\n', T_center_test_final, err_center(ii,jj), cpu_time(ii,jj));
    end
end

% Results table
fprintf('\nGrid and time-step study (error in T_center at t = %.2f s)\n', t_final);
    fprintf(' Nx = Ny |   dt   |  error [K] |  CPU time [s]\n');
fprintf('-------------------------------------------------\n');
for ii = 1:nNx
    for jj = 1:nDt
        fprintf(' %6d | %.3f | %11.3f | %12.3f\n', ...
            Nx_list(ii), dt_list(jj), 100*err_center(ii,jj), cpu_time(ii,jj));
    end
end

function [T_final, x, y, t_hist, T_center_hist, Q_bottom_hist, cpu_time] = HeatTransferAnalysis(Nx, Ny, dt, t_final, params, do_plots)

    % Anaylze parameters
    rho       = params.rho;
    cp        = params.cp;
    k         = params.k;
    Lx        = params.Lx;
    Ly        = params.Ly;
    T_initial = params.T_initial;
    T_bottom  = params.T_bottom;
    T_top     = params.T_top;
    T_left    = params.T_left;
    T_right   = params.T_right;
    alpha = k/(rho*cp);

    % Geometry & grid
    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    x  = linspace(0,Lx,Nx);
    y  = linspace(0,Ly,Ny);

    % Time
    Nt = round(t_final/dt);
    Fo_x = alpha*dt/dx^2; % fourier
    Fo_y = alpha*dt/dy^2; % fourier

    % Initial temperature field BCs at t = 0
    T = T_initial * ones(Ny, Nx);
    T(1,:)   = T_bottom;   % bottom (y = 0)
    T(Ny,:)  = T_top;      % top   (y = Ly)
    T(:,1)   = T_left;     % left  (x = 0)
    T(:,Nx)  = T_right;    % right (x = Lx)

    % CN matrices A and B
    Ntot = Nx * Ny;
    A = spalloc(Ntot, Ntot, 5*Ntot);
    B = spalloc(Ntot, Ntot, 5*Ntot);
    b = zeros(Ntot, 1);

    index = @(i,j) (i-1)*Ny + j;  % (i: x, j: y)

    for i = 1:Nx
        for j = 1:Ny
            p = index(i,j);

            % boundaries
            if i == 1 || i == Nx || j == 1 || j == Ny
                % Dirichlet node: T^{n+1} = T_bc
                A(p,p) = 1.0;

                if j == 1       % bottom
                    b(p) = T_bottom;
                elseif j == Ny  % top
                    b(p) = T_top;
                elseif i == 1   % left
                    b(p) = T_left;
                elseif i == Nx  % right
                    b(p) = T_right;
                end
                % B row stays zero (no T^n dependence in BC eq)

            else
                % Interior node: 2D Crank-Nicolson
                % A * T^{n+1} = B * T^n + b

                % A coefficients
                A(p,p)                 = 1 + Fo_x + Fo_y;
                A(p, index(i+1, j))    = -Fo_x/2;
                A(p, index(i-1, j))    = -Fo_x/2;
                A(p, index(i, j+1))    = -Fo_y/2;
                A(p, index(i, j-1))    = -Fo_y/2;

                % B coefficients
                B(p,p)                 = 1 - Fo_x - Fo_y;
                B(p, index(i+1, j))    =  Fo_x/2;
                B(p, index(i-1, j))    =  Fo_x/2;
                B(p, index(i, j+1))    =  Fo_y/2;
                B(p, index(i, j-1))    =  Fo_y/2;
            end
        end
    end

    % histories
    t_hist        = zeros(Nt,1); % Time past
    Q_bottom_hist = zeros(Nt,1); % Heat transfer
    T_center_hist = zeros(Nt,1); % Temperature

    i_center = ceil(Nx/2); % Center distribution
    j_center = ceil(Ny/2); % Center distribution

    T_vec = T(:);

    tic; % Elapsed time

    for n = 1:Nt
        t = n * dt;

        % RHS = B*T^n + b
        RHS   = B * T_vec + b;

        % solve A*T^{n+1} = RHS
        T_vec = A \ RHS;

        % re-impose Dirichlet BCs (robustness)
        % bottom
        for i = 1:Nx
            p_bottom = index(i,1);
            T_vec(p_bottom) = T_bottom;
        end
        % top
        for i = 1:Nx
            p_top = index(i,Ny);
            T_vec(p_top) = T_top;
        end
        % left / right
        for j = 1:Ny
            p_left  = index(1,j);
            p_right = index(Nx,j);
            T_vec(p_left)  = T_left;
            T_vec(p_right) = T_right;
        end

        T_mat = reshape(T_vec, Ny, Nx); % Temperature matrix for all grid

        % bottom heat transfer
        dTdy_bottom = (T_mat(2,:) - T_mat(1,:)) / dy;   % [K/m]
        q_bottom    = -k * dTdy_bottom;                 % [W/m^2]
        Qdot_bottom = sum(q_bottom) * dx;               % [W] per unit thickness

        % center temperature
        T_center = T_mat(j_center, i_center);

        % histories
        t_hist(n)        = t;
        Q_bottom_hist(n) = Qdot_bottom;
        T_center_hist(n) = T_center;
    end

    cpu_time = toc;

    % final field
    T_final = reshape(T_vec, Ny, Nx);

    % Eğer istersen base case için grafikler
    if do_plots
        
        [X, Y] = meshgrid(x, y);

        % final contour
        figure(1);
        contourf(X, Y, T_final, 50, 'LineColor', 'none');
        colorbar;
        colormap(gca, 'turbo');
        xlabel('x [m]');
        ylabel('y [m]');
        title(sprintf('Temperature distribution at t = %.2f s', t_final));
        axis equal tight;
        set(gcf, 'Color', 'w')
        caxis([300 400]);

        % center temperature vs time
        figure(2);
        plot(t_hist, T_center_hist, 'LineWidth', 1.5);
        xlabel('time [s]');
        ylabel('T_{center} [K]');
        title('Temperature at plate center vs. time');
        set(gcf, 'Color', 'w')
        grid on;

        % bottom heat transfer vs time
        figure(3);
        plot(t_hist, Q_bottom_hist, 'LineWidth', 1.5);
        xlabel('time [s]');
        ylabel('Q_{bottom} [W per unit thickness]');
        title('Total heat transfer rate through bottom boundary');
        set(gcf, 'Color', 'w')
        grid on;

        % mid-height profile
        j_mid = ceil(Ny/2);
        T_midline = T_final(j_mid,:);

        figure(4);
        plot(x, T_midline, 'LineWidth', 1.5);
        xlabel('x [m]');
        ylabel(sprintf('T(x, y = %.3f m) [K]', y(j_mid)));
        title('Temperature profile along x at mid-height (final time)');
        grid on;
        set(gcf, 'Color', 'w')
        drawnow
    end

end
