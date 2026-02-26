% MTH 261 — Coding Homework 2 

function d = detTri(U)
    [m, n] = size(U);
    if m ~= n
        error("Input is not a square matrix.")
    end

    tol = 1e-12;

    isUpper = true;
    isLower = true;
    
    for i = 1:n
        for j = 1:n
            if i > j && abs(U(i,j)) > tol 
                isUpper = false;
            end
            if i < j && abs(U(i,j)) > tol 
                isLower = false;
            end
        end
    end
    
    if ~isUpper && ~isLower
        error("Input matrix is not triangular.");
    end

    d = 1; % Initialize determinant
    for i = 1:n
        d = d * U(i, i); % Multiply diagonal entries
    end
end


function Uinv = invUpper(U)
    [m, n] = size(U);
    if m ~= n            % Check whether the matrix is a square
        error('Input matrix U is not a square.');
    end

    tol = 1e-12;

    isUpper = true;
    for i = 1:n
        for j = 1:n
            if i > j && abs(U(i,j)) > tol 
                isUpper = false;
            end
        end
    end

    if ~isUpper
        error("Input matrix is not upper triangular.");
    end

    if any(diag(U) == 0) % Check the diagonal entry
        error('Matrix U has a zero diagonal entry and the inverse does not exist.');
    end
    
    n = size(U, 1);
    Uinv = zeros(n);  % initialize an empty matrix that has the same size as U

    for i = 1:n
        ei = zeros(n, 1); 
        ei(i) = 1;          % a row for the identity matrix
        xi = backsub_local(U, ei); % solve U*x_i = e_i
        Uinv(:, i) = xi; % put the solution into the ith column
    end
end

% ---------------- OPTIONAL: local helper stub ----------------
% You may implement a small back-substitution helper inside this file
% (recommended). Do NOT use \ or built-in solvers.
%
% Tip: back-substitution solves from i=n down to 1.
%
function x = backsub_local(U,b)
    n = size(U, 1); % get the number of rows and columns in the matrix U
    x = zeros(n, 1); % initialize a row for the vector
    for i=n:-1:1 % the number of the row from the last one to the first one
        sum=0;   % initialize the sum of the known variables
        for j=i+1:n    % find the position of the known variables
            sum=sum+U(i,j)*x(j);  % for each row, add the sum of the known terms
        end
        x(i)=(b(i)-sum)/U(i,i);  % calculate the value of x with the sum of the known terms and the whole sum
    end
end


% ---------------- OPTIONAL: quick self-checks ----------------
% You may uncomment this block to test locally.
%
%
% % --- detTri checks ---
U1 = [2 1 -3; 0 4 5; 0 0 -1];
L1 = [3 0 0; -2 5 0; 1 7 2];
detU1 = detTri(U1)
detL1 = detTri(L1)
% % Compare (outside functions only):
det(U1), det(L1)
%
% % --- invUpper checks ---
U2 = [1 2; 0 2];
U3 = [2 -1 0; 0 3 4; 0 0 5];
U2inv = invUpper(U2)
U3inv = invUpper(U3)
norm(U2*U2inv - eye(2))
norm(U3*U3inv - eye(3))
% %
