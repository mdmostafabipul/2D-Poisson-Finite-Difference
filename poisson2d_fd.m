function [x, y, U, U_exact, err_inf, err_l2, A] = poisson2d_fd(N, L)
%POISSON2D_FD Solve -Delta u = f on (0,L)^2 using the 5-point stencil.
%
% Boundary condition:
%   u = 0 on the boundary
%
% Exact solution:
%   u(x,y) = sin(pi*x/L)*sin(pi*y/L)
%
% Input:
%   N  - total number of grid points in each coordinate direction
%   L  - side length of the square domain
%
% Outputs:
%   x, y      - full grid coordinates
%   U         - numerical solution on the full grid
%   U_exact   - exact solution on the full grid
%   err_inf   - maximum norm error
%   err_l2    - discrete L2 error
%   A         - sparse finite-difference matrix

if N < 4
    error('N must be at least 4.');
end

h = L/(N-1);
x = linspace(0, L, N);
y = linspace(0, L, N);

% Number of interior grid points in one direction
n = N-2;

%% Interior grid
xi = x(2:end-1);
yi = y(2:end-1);
[XI, YI] = meshgrid(xi, yi);

%% Right-hand side
% Exact: sin(pi*x/L) sin(pi*y/L)
% For -Delta u = f:
% f = 2*(pi/L)^2*sin(pi*x/L)*sin(pi*y/L)
F = 2*(pi/L)^2 .* sin(pi*XI/L) .* sin(pi*YI/L);
b = F(:);

%% 1D second-difference matrix for -d^2/dx^2
e = ones(n,1);
T = spdiags([-e, 2*e, -e], -1:1, n, n);

%% 2D five-point Laplacian matrix
I = speye(n);
A = (kron(I,T) + kron(T,I)) / h^2;

%% Solve the sparse linear system
u_vec = A\b;

%% Put the interior solution back on the full grid
U = zeros(N,N);
U(2:end-1,2:end-1) = reshape(u_vec, n, n);

%% Exact solution on the full grid
[X, Y] = meshgrid(x, y);
U_exact = sin(pi*X/L).*sin(pi*Y/L);

%% Errors
E = U-U_exact;
err_inf = max(abs(E(:)));
err_l2 = sqrt(h^2 * sum(E(:).^2));

end
