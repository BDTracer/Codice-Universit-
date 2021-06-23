function s3 = spline0(xi, yi, xq)
n = length(xi);
a = yi(1 : end - 1);         %Copiamo le ordinate nel vettore 
b = zeros(n - 1, 1);
d = zeros(n - 1, 1);

hi = diff(xi);
fhi = diff(yi);

A = zeros(n);
B = zeros(n, 1);
A(1, 1) = 1;
A(n, n) = 1;

% syms x;
% if (nargin == 2)
%     xx = sym("x",[1,n]);
% end

for i = 2 : n - 1
    A(i, i - 1) = hi(i - 1);
    A(i, i) = 2*(hi(i - 1) + hi(i));
    A(i, i + 1) = hi(i);
    B(i) = 3*(fhi(i) / hi(i) - fhi(i - 1) / hi(i - 1));
end

c = A \ B;

for i = 1 : n - 1
    d(i) = (c(i + 1) - c(i)) / (3 * hi(i));
    b(i) = fhi(i) / hi(i) - hi(i)*(2*c(i) + c(i + 1)) / 3;
end

[mxq, nxq] = size(xq);
s3 = zeros(mxq, nxq);

for i=1:mxq*nxq
    for k=1:n-1
        if xq(i) >= xi(k) && xq(i) < xi(k + 1)
            j = k;
            break;
        elseif xq(i) == xi(n)
            j = n-1;
        end
    end
    
    s3(i) = a(j) + b(j)*(xq(i) - xi(j)) + c(j)*(xq(i) - xi(j))^2 + d(j)*(xq(i) - xi(j))^3;
end

