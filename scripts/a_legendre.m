function [P] = a_legendre(a, b, x)
%a_legendre Calculates associated Legendre Polynomial value
%   Detailed explanation goes here

syms y;
f = power(y.^2 - 1, a);
for i = 1:(a + b)
    f = diff(f);
end

g = (-1)^b * power(1 - y.^2, b/2) / (2^a * factorial(a));

P = vpa(subs(g.*f, y, x));

end
