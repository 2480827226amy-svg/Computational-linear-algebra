% MTH 261 — Coding Homework 1 (Template)
% Groupmate: Amy Tang, Selina Fang, Belle Wen


function x = backsub(U,b)
%BACKSUB  Solve Ux=b when U is upper triangular (no backslash).
%
% Inputs:
%   U : (n x n) upper triangular matrix (assume nonzero diagonal)
%   b : (n x 1) column vector
%
% Output:
%   x : (n x 1) solution vector

    % Hint (setup): n = size(U,1);  x = zeros(n,1);

    % Common bug: loop direction matters (back-substitution goes from n down to 1).

    % --- TODO: back-substitution loop ---
    % Implement the algorithm to solve Ux=b without using MATLAB solvers.
    % Suggested structure (do not copy/paste blindly; make it yours):
    %   - loop i from n down to 1
    %   - compute the running sum of known terms
    %   - solve for x(i)

    %the number of rows and columns in the matrix U
    n = size(U, 1); 

    %the exact row we want in the matrix
    x = zeros(n, 1); 

    %the number of the row from the last one to the first one
    for i=n:-1:1 

        %initialize the sum of the known variables
        sum=0; 

        %the position of the known variables
        for j=i+1:n 

            %for each row, add the sum of the known terms
            sum=sum+U(i,j)*x(j); 
        end

        %calculate the value of x with the sum of the know terms and the whole sum
        x(i)=(b(i)-sum)/U(i,i); 
    end
end


function A = FirstPivot(A)
%FIRSTPIVOT  First pivot step of Gaussian elimination (no pivoting).
%
% Input:
%   A : (n x (n+1)) augmented matrix (assume A(1,1) ~= 0)
%
% Output:
%   A : matrix after:
%       1) scaling row 1 so A(1,1)=1
%       2) clearing below pivot in column 1

    % Hint (setup): [n,~] = size(A);  p = A(1,1);

    % --- TODO: scale first row ---
    % Scale row 1 so that A(1,1) becomes 1.

    % --- TODO: eliminate below pivot ---
    % For each row i=2,...,n, eliminate A(i,1) using a row operation.

   [n, ~] = size(A);
    % Scale row 1 so that A(1,1) becomes 1
    p = A(1,1);

    % Eliminate entries below the pivot in column 1
    A(1,:) = A(1,:)/p; 

    for i = 2:n
        A(i,:) = A(i,:) - A(i,1)*A(1,:);
    end
end


% ---------------- OPTIONAL: quick self-checks ----------------
% You may UNCOMMENT this section by removing comments between %{ and %}
% to test locally, but do not rely on it
% for grading. Grading will call the functions above on several inputs.
%
%
% % Test backsub:
U = [1 5; 0 2];
b = [4; 1];
x = backsub(U,b)
norm(U*x - b)  

% % Test FirstPivot:
A = [ 2  1 -1   8;
     -2 -1  2 -11;
     -2  1  2  -3];
Aout = FirstPivot(A)
Aout(:,1)