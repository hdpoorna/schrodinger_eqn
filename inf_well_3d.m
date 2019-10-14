% inf_well_3d.m : infinitely deep potential well in 3D
% author : hdpoorna
% MATLAB R2018b

%% Initialization

clc;
clear;
close all;

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);

m_e = 9.109e-31;    % mass of an electron
m = m_e;

Lx = input('Enter Lx (>0): ');
Ly = input('Enter Ly (>0): ');
Lz = input('Enter Lz (>0): ');

x = 0:Lx/100:Lx;
y = 0:Ly/100:Ly;
z = 0:Lz/100:Lz;
[X,Y,Z] = meshgrid(x,y,z);

nx = input('Enter nx (1,2,..): ');
ny = input('Enter ny (1,2,..): ');
nz = input('Enter nz (1,2,..): ');

%% Plotting probability distribution

figure('Name', 'Probability Distribution');
psi_xyz = sqrt(8/(Lx*Ly*Lz))*sin(nx*pi*X/Lx).*sin(ny*pi*Y/Ly).*sin(nz*pi*Z/Lz);

prob_xyz = abs(psi_xyz).^2;

% slice 1
xslice = Lx/2;        % location of y-z planes
yslice = Ly/2;        % location of x-z planes
zslice = Lz/2;        % location of x-y planes

%{
% slice 3
xslice = [Lx/4, Lx/2, 3*Lx/4];        % location of y-z planes
yslice = [Ly/4, Ly/2, 3*Ly/4];        % location of x-z planes
zslice = [Lz/4, Lz/2, 3*Lz/4];        % location of x-y planes
%}

slice(X, Y, Z, prob_xyz, xslice, yslice, zslice)
xlabel('x'); ylabel('y'); zlabel('z'); colorbar; daspect([1 1 1]);
title(sprintf('Probability Distribution for n_x = %s, n_y = %s, n_z = %s', num2str(nx), num2str(ny), num2str(nz)));