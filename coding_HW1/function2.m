function A = FirstPivot(A)
    [n, ~] = size(A);

    p = A(1,1);
    A(1,:) = A(1,:)/p;

    for i = 2:n; j
        A(i,:) = A(i,:) - A(i,1)*A(1,:);
    end
end


A = [ 2  1 -1   8;
     -3 -1  2 -11;
     -2  1  2  -3];
Aout = FirstPivot(A)
Aout(:,1)