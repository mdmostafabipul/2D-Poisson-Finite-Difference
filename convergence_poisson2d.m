%% Grid Convergence Study for the 2D Poisson Equation
clear; clc; close all;

L = 1.0;

% Total grid points per direction
N_values = [11, 21, 41, 81, 161];

h_values = zeros(size(N_values));
errors_inf = zeros(size(N_values));
errors_l2  = zeros(size(N_values));
unknowns   = zeros(size(N_values));
solve_time = zeros(size(N_values));

for k = 1:length(N_values)

    N = N_values(k);
    h_values(k) = L/(N-1);
    unknowns(k) = (N-2)^2;

    num_runs = 10;
times = zeros(num_runs,1);

for j = 1:num_runs

    tic;

    [~, ~, ~, ~, err_inf, err_l2] = ...
        poisson2d_fd(N, L);

    times(j) = toc;

end

solve_time(k) = mean(times);

    errors_inf(k) = err_inf;
    errors_l2(k)  = err_l2;

end

%% Print convergence table
fprintf(' N        h          Unknowns      Linf Error       Order\n');
fprintf('-------------------------------------------------------------\n');

for k = 1:length(N_values)

    if k == 1
        fprintf('%3d   %.4e   %8d    %.4e        ---\n', ...
            N_values(k), h_values(k), unknowns(k), errors_inf(k));
    else
        p = log(errors_inf(k-1)/errors_inf(k)) / ...
            log(h_values(k-1)/h_values(k));

        fprintf('%3d   %.4e   %8d    %.4e      %.4f\n', ...
            N_values(k), h_values(k), unknowns(k), ...
            errors_inf(k), p);
    end

end

fprintf('\nAverage Runtime Results\n');
fprintf('------------------------------\n');

for k = 1:length(N_values)

    fprintf('N = %3d, Unknowns = %6d, Average time = %.6e s\n', ...
        N_values(k), unknowns(k), solve_time(k));

end

%% Error convergence plot
figure;
loglog(h_values, errors_inf, 'o-', 'LineWidth', 1.8);
xlabel('h');
ylabel('L_\infty error');
title('Grid Convergence for the 2D Poisson Equation');
grid on;
set(gca,'XDir','reverse');

%% Runtime plot
figure;
loglog(unknowns, solve_time, 'o-', 'LineWidth', 1.8);
xlabel('Number of unknowns');
ylabel('Solve time (seconds)');
title('Runtime versus Problem Size');
grid on;
