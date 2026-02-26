% MTH 261 — Coding Homework 3
% Groupmate: Amy Tang, Selina Fang, Belle Wen

%% ------------------- Script section (required) -------------------
% Solve Ax=b using YOUR LU factorization and substitution functions.

% Given system:
A = [  7  -2   1;
      14  -7  -3;
      -7  11  18];
b = [12; 17; 5];

% --- TODO: compute LU factorization and solve ---
[L,U] = LU(A);
y = fwdsub(L,b);
x = backsub(U,y);

% Display your result x:
disp(x);

% Recommended numerical checks (OK to use built-ins here):
fprintf('||L*U - A|| = %.6g\n', norm(L*U - A));
fprintf('||A*x - b|| = %.6g\n', norm(A*x - b));
x_check = A\b;  % allowed for checking only
fprintf('||x - x_check|| = %.6g\n', norm(x - x_check));


%% ------------------- Function stubs (to implement) -------------------

function [L, U] = LU(A)
    n = size(A, 1); % get matrix dimension
    L = eye(n); % initialize L as identity (unit lower triangular)
    U = A; % initialize U as a copy of A
    
    for k = 1:n-1 % iterate over each pivot column

        for i = k+1:n % iterate over rows below the pivot

            % Compute multiplier for row i
            m = U(i,k) / U(k,k);
            L(i,k) = m;

            % Eliminate column k in rows below pivot
            U(i,:) = U(i,:) - m * U(k,:);
        end
    end
end


function y = fwdsub(L,b)
    n = size(L, 1);
    y = zeros(n, 1);  % initialize solution vector y with zeros
    for i = 1:n    % loop from first row to last row 
        y(i) = b(i) - L(i, 1:i-1) * y(1:i-1)/ L(i, i);
    end
end


function x = backsub(U, y)
    % Solve Ux = y using back substitution.
    % U is upper triangular, so we solve from the last row upward.

    n = length(y);          % Number of unknowns
    x = zeros(n, 1);        % Initialize solution vector
    
    for i = n:-1:1
        % Compute their contribution and solve for x(i).
        x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
    end
end


% ---------------- OPTIONAL: quick self-checks ----------------
% You may uncomment this block to test pieces locally.
%
% % Quick test for fwdsub on a simple L:
Ltest = [1 0 0; 2 1 0; -1 4 1];
btest = [1; 0; 3];
ytest = fwdsub(Ltest,btest)
norm(Ltest*ytest - btest)
%
% % Quick test for backsub on a simple U:
Utest = [1 2 3; 0 2 -1; 0 0 5];
btest2 = [1; 4; 10];
xtest = backsub(Utest,btest2)
norm(Utest*xtest - btest2)

