% H_atom.m - Hydrogen atom orbitals
% author : hdpoorna
% MATLAB R2018b

%% Initialization

clc;
clear;
close all;

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);
eps = 8.854e-12;    % Permittivity of free space

q = 1.602e-19;      % charge of an electron
m_p = 1.6726e-27;   % mass of a proton
m_e = 9.109e-31;    % mass of an electron
mass = m_e;         % assuming proton is fixed, neglecting reduced mass. (0.9995m_e)
a0 = 4*pi*eps*power(h_bar/q, 2)/mass;   % Bohr radius

n = input('Enter n (>= 1 ; a positive integer): ');      % >= 1
l = input('Enter l (0 <= l <= n-1 ; a non-negative integer): ');      % 0 <= l < n
m = input('Enter m (-l <= m <= +l ; an integer): ');      % -l <= m <= l

%% Plotting energy levels

n_max = 5;
z_E = 0:0.01:1;

E = zeros(n_max, length(z_E));

figure('Name', 'Energy Levels'); hold on;
for nE = 1:n_max
    E_nE = -mass*power((q^2)/(pi*eps*h_bar*nE), 2)/32;
    E(nE, :) = E_nE;
    plot(z_E, E(nE, :), 'r-');
    text(max(z_E), max(E(nE, :)), sprintf('n = %s', num2str(nE)));
end
hold off; grid on;
title('Energy Levels'); ylabel('Energy');

%% Plotting probability distributions and orbital shapes

r_max = 120e-12;        % Van der Vaal's radius of Hydrogen
dp = 15;

x = [-r_max:r_max/dp:-r_max/dp r_max/dp:r_max/dp:r_max];
y = [-r_max:r_max/dp:-r_max/dp r_max/dp:r_max/dp:r_max];
z = [-r_max:r_max/dp:-r_max/dp r_max/dp:r_max/dp:r_max];

[X, Y, Z] = meshgrid(x, y, z);

R = sqrt(X.^2 + Y.^2 + Z.^2);
Theta = atan2(sqrt(X.^2 + Y.^2), Z);
Phi = atan2(Y,X);

%{
Rnl_r = sqrt(factorial(n - l - 1)/(2*n*(factorial(n + 1))^3))*power(2/(n*a0), l + 3/2)*...
    a_laguerre(n + l, 2*l + 1, 2*R/(n*a0)).*...
    power(R, l).*exp(-R/(n*a0));

Ylm = sqrt((2*l + 1)*factorial(l - abs(m))/(4*pi*factorial(l + abs(m))))*...
    a_legendre(l, abs(m), cos(Theta)).*...
    exp(1i*m*Phi);
%}

psi = H_wave_f(n, l, m, R, Theta, Phi);
prob = double(abs(psi).^2);

% slice 1
xslice = 0;        % location of y-z planes
yslice = 0;        % location of x-z planes
zslice = 0;        % location of x-y planes

%{
% slice 3
xslice = [-r_max/4, 0, r_max/4];        % location of y-z planes
yslice = [-r_max/4, 0, r_max/4];        % location of x-z planes
zslice = [-r_max/4, 0, r_max/4];        % location of x-y planes
%}


figure('Name', 'Probability Distribution');
subplot(1,2,1);
slice(X, Y, Z, prob, xslice, yslice, zslice)
xlabel('x'); ylabel('y'); zlabel('z'); colorbar; daspect([1 1 1]);
title(sprintf('Probability Distribution for n = %s, l = %s, m = %s', num2str(n), num2str(l), num2str(m)));


%figure('Name', 'Orbital Shape');
subplot(1,2,2);
iso_val = max(prob, [], 'all') - (max(prob, [], 'all') - min(prob, [], 'all'))*0.85;
p = patch(isosurface(X, Y, Z, prob, iso_val));
p.FaceColor = 'red'; p.EdgeColor = 'None';
daspect([1 1 1]); view(3); camlight;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Orbital Shape for n = %s, l = %s, m = %s', num2str(n), num2str(l), num2str(m)));

%% Define local functions

function [psi] = H_wave_f(n, l, m, R, Theta, Phi)
%h_wave_f Calculates Hydrogen atom wave function

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


function [L] = a_laguerre(a, b, x)
%a_laguerre Calculates associated Laguerre Polynomial value

syms y;
f = power(y, a + b).*exp(-y);
for i = 1:a
    f = diff(f);
end

g = power(y, -b).*exp(y)/factorial(a);

L = vpa(subs(g.*f, y, x));

end


function [P] = a_legendre(a, b, x)
%a_legendre Calculates associated Legendre Polynomial value

syms y;
f = power(y.^2 - 1, a);
for i = 1:(a + b)
    f = diff(f);
end

g = (-1)^b * power(1 - y.^2, b/2) / (2^a * factorial(a));

P = vpa(subs(g.*f, y, x));

end
