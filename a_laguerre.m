function [L] = a_laguerre(a, b, x)
%a_laguerre Calculates associated Laguerre Polynomial value
%   Detailed explanation goes here

syms y;
f = power(y, a + b).*exp(-y);
for i = 1:a
    f = diff(f);
end

g = power(y, -b).*exp(y)/factorial(a);

L = vpa(subs(g.*f, y, x));

end
