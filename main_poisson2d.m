%% 2D Poisson Equation using a 5-point finite difference method
% Problem:
%   -Delta u = f,  (x,y) in (0,1)x(0,1)
%   u = 0 on the boundary
%
% Exact solution:
%   u(x,y) = sin(pi*x)*sin(pi*y)
%
% Therefore:
%   f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)

clear; clc; close all;

%% Parameters
N = 41;                 % total grid points in each direction
L = 1.0;                % domain length

%% Solve the problem
[x, y, U_num, U_ex, err_inf, err_l2, A] = poisson2d_fd(N, L);

%% Display basic information
h = L/(N-1);
nInterior = N-2;

fprintf('2D Poisson Equation using 5-point finite differences\n');
fprintf('---------------------------------------------------\n');
fprintf('Total grid points per direction N = %d\n', N);
fprintf('Interior points per direction    = %d\n', nInterior);
fprintf('Grid spacing h                   = %.6e\n', h);
fprintf('Number of unknowns               = %d\n', nInterior^2);
fprintf('Nonzeros in sparse matrix A      = %d\n', nnz(A));
fprintf('\nErrors:\n');
fprintf('L-infinity error                 = %.6e\n', err_inf);
fprintf('Discrete L2 error                = %.6e\n', err_l2);

%% Plot numerical solution
[X, Y] = meshgrid(x, y);

figure;
surf(X, Y, U_num);
xlabel('x');
ylabel('y');
zlabel('u_h(x,y)');
title('Numerical Solution of the 2D Poisson Equation');
shading interp;
view(45,30);

%% Plot exact solution
figure;
surf(X, Y, U_ex);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Exact Solution');
shading interp;
view(45,30);

%% Plot absolute error
figure;
surf(X, Y, abs(U_num-U_ex));
xlabel('x');
ylabel('y');
zlabel('|u_h-u|');
title('Pointwise Absolute Error');
shading interp;
view(45,30);

%% Compare through the line y = 0.5
[~, mid_idx] = min(abs(y-0.5));

figure;
plot(x, U_ex(mid_idx,:), 'LineWidth', 1.8);
hold on;
plot(x, U_num(mid_idx,:), 'o', 'MarkerSize', 4);
xlabel('x');
ylabel('u(x,0.5)');
title('Cross-section at y = 0.5');
legend('Exact solution', 'Numerical solution', 'Location', 'best');
grid on;

figure;
spy(A);
title('Sparsity Pattern of the 2D Poisson Matrix');
xlabel('Column Index');
ylabel('Row Index');