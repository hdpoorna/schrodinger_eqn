% inf_well_1d.m : Finitely deep potential well in 1D
% author : hdpoorna
% MATLAB R2018b

%% Initialization

clc;
clear;
close all;

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);
Lz = input('Enter Lz (>0): ');
m_e = 9.109e-31;    % mass of an electron
m = m_e;

% E1_inf = (h_bar*pi/Lz)^2 / (2*m);
% V0 = 10*E1_inf;
v0 = input('Enter v0 (>=1.5) : ');      % = V0/E1_inf
% for better visualization v0 >= 1.5

eps = linspace(0, v0, 100*v0 + 1);

z1 = -3*Lz/2:Lz/100:-Lz/2 -Lz/100;
z2 = -Lz/2 + Lz/100:Lz/100:Lz/2 - Lz/100;
z3 = Lz/2 + Lz/100:Lz/100:3*Lz/2;
z = [z1 z2 z3];

%% for symmetric solution

y1_sym = sqrt(eps).*tan(sqrt(eps)*pi/2);
y2_sym = sqrt(v0 - eps);

inds_sym = find(mod(sqrt(eps), 2) == 1);        % find discontinuities
y1_sym(inds_sym) = NaN;
y2_sym(inds_sym) = NaN;
eps(inds_sym) = NaN;

figure('Name', 'Finding Energy Levels for Symmetric Solution');
plot(eps, y1_sym, 'r-', eps, y2_sym, 'b-');
title(sprintf('Finding Energy Levels for Symmetric Solution, when v_0 = %s', num2str(v0)));
xlabel('\epsilon'); ylabel('y'); ylim([0 sqrt(v0)]); hold on; grid on; 
[eps_i, yi] = polyxpoly(eps, y1_sym, eps, y2_sym);      % finding intersections
scatter(eps_i, yi, 'ko');
legend('LHS', 'RHS');


cL = cos(pi*sqrt(eps_i)/2)./exp(-pi*sqrt(v0 - eps_i)/2);
B = sqrt(2./(1 +...
    cL.^2 .* (exp(pi*sqrt(v0 - eps_i)) + exp(-pi*sqrt(v0 - eps_i)))./(pi*sqrt(v0 - eps_i)) +...
    sin(pi*sqrt(eps_i))./(pi*sqrt(eps_i))));


psi1_sym = B.*cL.*exp(pi*sqrt(v0 - eps_i).*z1/Lz);
psi2_sym = B.*cos(pi*sqrt(eps_i).*z2/Lz);
psi3_sym = B.*cL.*exp(-pi*sqrt(v0 - eps_i).*z3/Lz);
psi_sym = [psi1_sym psi2_sym psi3_sym];


figure('Name', 'Wave Function and Probability for Symmetric Solution'); hold on;
for ni = 1:length(eps_i)
    subplot(length(eps_i), 2, 2*length(eps_i) + 1 - 2*ni);
    plot(z, psi_sym(ni, :), 'r-'); grid on;
    title(sprintf('Wave function for \\epsilon = %s when v_0 = %s', num2str(eps_i(ni)), num2str(v0)));
    xlabel('z'); ylabel(sprintf('\\psi_{%s}(z)', num2str(ni)));
    
    subplot(length(eps_i), 2, 2*length(eps_i) + 2 - 2*ni);
    plot(z, abs(psi_sym(ni, :)).^2, 'b-'); grid on;
    title(sprintf('Probability Distribution for \\epsilon = %s when v_0 = %s', num2str(eps_i(ni)), num2str(v0)));
    xlabel('z'); ylabel(sprintf('|\\psi_{%s}(z)|^{2}', num2str(ni)));
end

%% For antisymmetric solution


y1_asym = -1*sqrt(eps).*cot(sqrt(eps)*pi/2);
y2_asym = sqrt(v0 - eps);

inds_asym = find(mod(sqrt(eps), 2) == 0);       % finding discontinuities
y1_asym(inds_asym) = NaN;
y2_asym(inds_asym) = NaN;
eps(inds_asym) = NaN;

figure('Name', 'Finding Energy Levels for Antisymmetric Solution');
plot(eps, y1_asym, 'r-', eps, y2_asym, 'b-');
title(sprintf('Finding Energy Levels for Antisymmetric Solution, when v_0 = %s', num2str(v0)));
xlabel('\epsilon'); ylabel('y'); ylim([0 sqrt(v0)]); hold on; grid on;
[eps_j, yj] = polyxpoly(eps, y1_asym, eps, y2_asym);        % finding intersections
scatter(eps_j, yj, 'ko');
legend('LHS', 'RHS');

sL = sin(pi*sqrt(eps_i)/2)./exp(-pi*sqrt(v0 - eps_i)/2);
A = sqrt(2./(1 +...
    sL.^2 .* (exp(pi*sqrt(v0 - eps_i)) + exp(-pi*sqrt(v0 - eps_i)))./(pi*sqrt(v0 - eps_i)) -...
    sin(pi*sqrt(eps_i))./(pi*sqrt(eps_i))));


psi1_asym = -A.*sL.*exp(pi*sqrt(v0 - eps_i).*z1/Lz);
psi2_asym = A.*sin(pi*sqrt(eps_i).*z2/Lz);
psi3_asym = A.*sL.*exp(-pi*sqrt(v0 - eps_i).*z3/Lz);
psi_asym = [psi1_asym psi2_asym psi3_asym];


figure('Name', 'Wave Function and Probability for Antisymmetric Solution'); hold on;
for nj = 1:length(eps_j)
    subplot(length(eps_j), 2, 2*length(eps_j) + 1 - 2*nj);
    plot(z, psi_asym(nj, :), 'r-'); grid on;
    title(sprintf('Wave function for \\epsilon = %s when v_0 = %s', num2str(eps_j(nj)), num2str(v0)));
    xlabel('z'); ylabel(sprintf('\\psi_{%s}(z)', num2str(nj)));
    
    subplot(length(eps_j), 2, 2*length(eps_j) + 2 - 2*nj);
    plot(z, abs(psi_asym(nj, :)).^2, 'b-'); grid on;
    title(sprintf('Probability Distribution for \\epsilon = %s when v_0 = %s', num2str(eps_j(nj)), num2str(v0)));
    xlabel('z'); ylabel(sprintf('|\\psi_{%s}(z)|^{2}', num2str(nj)));
end
