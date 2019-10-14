% eigen_H - solve for eigen energy values and wavefunctions for H atom
% author : hdpoorna
% MATLAB R2018b

%% Initialization

clc;
clear;
close all;

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);
L = 1;
r = 0.5;
m_e = 9.109e-31;    % mass of an electron
eps = 8.854e-12;    % Permittivity of free space
q = 1.602e-19;      % charge of an electron

% k = q^2/(4*pi*eps);

c = (h_bar^2)/(2*m_e);
d = 1;

eigenV = input('For which eigen energy level do you want to plot? (1,2,..5) : ');

model = createpde;

%% Potential function

a = @(location,state)(-1./sqrt((location.x).^2 + (location.y).^2 + (location.z).^2));
% used potential function as (-1/r) for easier numerical calculations
% for MATLAB

%% Boundary condition

gm = multisphere(r);
model.Geometry = gm;
figure('Name', 'Boundary Conditions');
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5);
title('Boundary Conditions');
applyBoundaryCondition(model, 'dirichlet', 'Face', 1, 'u', 0);

%% Specify coefficients

specifyCoefficients(model, 'c', 1, 'a', a, 'd', d, 'm', 0, 'f', 0);
% c is normalized to 1 for easier numerical calculations for MATLAB

%% Solve PDE for eigenvalues

% generateMesh(model, 'Hmax', 0.05);
generateMesh(model);
evr = [-inf,200];
results = solvepdeeig(model, evr);

%% Plot probability distribution and orbital shape

[X,Y,Z] = meshgrid(-r:r/10:r, -r:r/10:r, -r:r/10:r);
V = interpolateSolution(results,X,Y,Z,eigenV);
V = reshape(V,size(X));
prob = abs(V).^2;

xslice = 0;        % location of y-z planes
yslice = 0;        % location of x-z planes
zslice = 0;        % location of x-y planes
figure('Name', 'Probability Distribution and Orbital Shape');
subplot(1,2,1);
slice(X, Y, Z, prob, xslice, yslice, zslice);
title(sprintf('Probability Distribution for eigen energy level %s', num2str(eigenV)));
xlabel('x'); ylabel('y'); zlabel('z'); colorbar; daspect([1 1 1]);

%figure('Name', 'Orbital Shape');
subplot(1,2,2);
iso_val = max(prob, [], 'all') - (max(prob, [], 'all') - min(prob, [], 'all'))*0.85;
p = patch(isosurface(X, Y, Z, prob, iso_val));
p.FaceColor = 'red'; p.EdgeColor = 'None';
title(sprintf('Orbital Shape for eigen energy level %s', num2str(eigenV)));
daspect([1 1 1]); view(3); camlight;
xlabel('x'); ylabel('y'); zlabel('z');