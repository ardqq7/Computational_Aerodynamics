% UCK419E - Computational Aerodynamics - Homework 2 - Question 3 Code
% Arda Ceker 110200131 
% Written in MATLAB R2025a

clear, clc, clf, close all, format long g;

data = readmatrix('input.txt'); % Reads file and saves data as a matrix
x = data(:,1); % Stores x values from first column
y = data(:,2); % Stores y values from first column

position = [0, 0]; % Starting from the origin
figure(1)
    hold on;
    for i = 1 : length(x)
        final = position + [x(i), y(i)]; % Adds start point of the element to others end point
        quiver(position(1), position(2), x(i), y(i), 0); % Draws vector line
        position = final; % Sets the new start point is the previous end point
    end

    % Plotting settings
    axis equal; grid minor;
    xlabel('X'); ylabel('Y');
    title('Sum of Vectors');
        set(gcf, 'Position', [100, 100, 1500, 800]); %, 'Color', 'w');
        set(findall(gca, 'Type', 'Text'), 'FontWeight', 'bold', 'FontName', 'Courier');
        set(findall(gca, 'Type', 'Quiver'), 'LineWidth', 2, 'MaxHeadSize', 0.2);
        set(gca, 'FontWeight', 'bold', 'FontName', 'Courier','FontSize', 14, 'LineWidth', 2);

[r, theta] = cartesian2polar(x, y);  % Calls the function cartesian2polar

output = [r, theta];
writematrix(output, 'output.txt', 'Delimiter', 'tab'); % Writes the output.txt file as a two column matrix

disp('Inout file is: (x , y)');
disp(data);
disp('Output file is: (r , θ)');
disp(output);

function [r, theta] = cartesian2polar(x, y)
    r = sqrt(x.^2 + y.^2); % r component 
    theta = atan2d(y, x); % theta component
end
