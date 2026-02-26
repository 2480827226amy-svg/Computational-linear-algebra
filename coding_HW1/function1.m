function x = backsub(U,b)
    n = size(U, 1);
    x = zeros(n, 1);
    for i=n:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+U(i,j)*x(j);
        end
    x(i)=(b(i)-sum)/U(i,i);
    end
end

U = [1 5; 2 2];
b = [4; 1];
x = backsub(U,b)
norm(U*x - b)