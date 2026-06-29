% UCK419E - Computational Aerodynamics - Homework 1 - Question 3 Code
% Arda Ceker 110200131 
% Written in MATLAB R2025a

clear, clc;
%  Input 
lb = input('Enter lower limit: ');
ub = input('Enter upper limit: ');

P = primeNumberFinder(lb, ub);

% Output 
disp('Prime numbers between the interval are:');
disp(P);

function P = primeNumberFinder(lb, ub)

    lb = ceil(lb); 
    ub = floor(ub);

    P = []; 

    for n = lb:ub
        if n < 2
            continue; % Skips for 1
        end

        test = true; % Test parameter. Returns false if the number is divided by other numbers.

        for d = 2:ub
            if d == n
                continue; % Skips to divide itself
            end

            if mod(n, d) == 0 % Division test
                test = false;
                break;
            end
        end

        if test
            P(end+1) = n;  % Adds into the matrix if it is prime
        end
    end
end
