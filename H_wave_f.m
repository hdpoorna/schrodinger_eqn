function [psi] = H_wave_f(n, l, m, R, Theta, Phi)
%h_wave_f Calculates Hydrogen atom wave function
%   Detailed explanation goes here

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);
eps = 8.854e-12;

q = 1.602e-19;
m_p = 1.2726e-27;
m_e = 9.109e-31;
mass = m_e;            % assuming proton is fixed, neglecting reduced mass. (o.9995m_e)
a0 = 4*pi*eps*power(h_bar/q, 2)/mass;

%{
Rnl_r = sqrt(factorial(n - l - 1)/(2*n*(factorial(n + 1))^3))*power(2/(n*a0), l + 3/2)*...
    a_laguerre(n + l, 2*l + 1, 2*R/(n*a0)).*...
    power(R, l).*exp(-R/(n*a0));
%}
Rnl_r = sqrt(factorial(n - l - 1)/(2*n*(factorial(n + 1))^3))*power(2/(n*a0), l + 3/2)*...
    a_laguerre(n + l, 2*l + 1, 2*R/(n*a0)).*...
    power(R, l).*exp(-R/(n*a0));

if m == 0
    Ylm = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
        a_legendre(l, abs(m), cos(Theta)).*...
        exp(1i*m*Phi);
    psi = Rnl_r.*Ylm;
elseif m > 0
    Ylm_pos = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
        a_legendre(l, abs(m), cos(Theta)).*...
        exp(1i*abs(m)*Phi);
    Ylm_neg = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
        a_legendre(l, abs(m), cos(Theta)).*...
        exp(-1i*abs(m)*Phi);
    psi = Rnl_r.*(Ylm_pos + Ylm_neg)/sqrt(2);
else
    Ylm_pos = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
        a_legendre(l, abs(m), cos(Theta)).*...
        exp(1i*abs(m)*Phi);
    Ylm_neg = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
        a_legendre(l, abs(m), cos(Theta)).*...
        exp(-1i*abs(m)*Phi);
    psi = -1i*Rnl_r.*(Ylm_pos - Ylm_neg)/sqrt(2);
end

end
