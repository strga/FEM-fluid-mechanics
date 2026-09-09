function [x,Y] = FEM_1D(n1, f)

%f = @(x) x + 2;

% Interval division

n = n1 - 1;             % Discretize to n values / nodes ... dim V = n
h = 1. / (n1);          % step size ... steps to get from v(0) to v(1) when I seek n-values in between v(0) and v(1) (= dim V)
x = (0:h:1)';           % discretization of x
y = f(x);               % final values from function f = @(x) 

% Stiffness matrix

U = diag(2 * ones(n,1));
V = diag(-ones(n-1,1),1);
W = diag(-ones(n-1,1),-1);
A = U + V + W;                  % Sparse stiffness matrix

b = h * h * f(x(2:end - 1));    % Load vector computation of inner nodes
                                % h = 1/(1+n) 1/h Au = b and also b = h*f(x)'
u = A \ b;

Y = [0;u;0];                    % Column vector with 1st and last value of 0 + inner nodes - BC enforcement

end

