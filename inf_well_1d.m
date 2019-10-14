% inf_well_1d.m : infinitely deep potential well in 1D
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

z = 0:Lz/1000:Lz;
n = 1:5;

E = zeros(length(n), length(z));
psi_z = zeros(length(n), length(z));
prob_z = zeros(length(n), length(z));

%% Plotting energy levels

figure('Name', 'Energy Levels'); hold on;
for nE = 1:length(n)
    E_nE = power((h_bar*nE*pi/Lz), 2)/(2*m);
    E(nE, :) = E_nE;
    plot(z, E(nE, :), 'r-');
    text(max(z), max(E(nE, :)), sprintf('n = %s', num2str(nE)));
end
hold off; grid on;
title('Energy Levels'); xlabel('z'); ylabel('Energy');

%% Plotting wave function and probability distributions

figure('Name', 'Wave Function and Probability'); hold on;
for nz = 1:length(n)
    psi_nz = sqrt(2/Lz)*sin(nz*pi*z/Lz);
    psi_z(nz, :) = psi_nz;
    subplot(length(n), 2, 2*length(n) + 1 - 2*nz);
    plot(z, psi_nz, 'b-'); grid on;
    title(sprintf('Wave function for n = %s', num2str(nz)));
    xlabel('z'); ylabel(sprintf('\\psi_{%s}(z)', num2str(nz)));
    text(max(z), max(psi_nz), sprintf('n = %s', num2str(nz)));
    
    prob_nz = abs(psi_nz).^2;
    prob_z(nz, :) = prob_nz;
    subplot(length(n), 2, 2*length(n) + 2 - 2*nz);
    plot(z, prob_nz, 'g-'); grid on;
    title(sprintf('Probability Distribution for n = %s', num2str(nz)));
    xlabel('z'); ylabel(sprintf('|\\psi_{%s}(z)|^{2}', num2str(nz)));
    text(max(z), max(prob_nz), sprintf('n = %s', num2str(nz)));
end
hold off;