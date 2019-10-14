% eigen_box - solve for eigen energy values and wavefunctions for
%               infinitely deep potential well in 3D
% author : hdpoorna
% MATLAB R2018b

%% Initialization

clc;
clear;
close all;

h = 6.626e-34;      % Planck's constant
h_bar = h/(2*pi);
L = 1;
m_e = 9.109e-31;    % mass of an electron

% c = (h_bar^2)/(2*m_e);
d = 1;

eigenV = input('For which eigen energy level do you want to plot? (1,2,..5) : ');

model = createpde;

a = 0;      % potential

%% Boundary conditions

[x,y,z] = meshgrid([-L/2 L/2]);
x = x(:);
y = y(:);
z = z(:);
K = convhull(x,y,z);
nodes = [x';y';z'];
elements = K';
geometryFromMesh(model,nodes,elements);
figure('Name', 'Boundary Conditions');
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5);
title('Boundary Conditions');
applyBoundaryCondition(model, 'dirichlet', 'Face', 1:6, 'u', 0);

%% Specify coefficients

specifyCoefficients(model, 'c', 1, 'a', a, 'd', d, 'm', 0, 'f', 0);
% c is normalized to 1 for easier numerical calculations for MATLAB

%% Solve PDE for eigenvalues

% generateMesh(model, 'Hmax', 0.05);
generateMesh(model);
evr = [1,200];
results = solvepdeeig(model, evr);

%{
pdeplot3D(model,'ColorMapData',results.NodalSolution(:,3),...
    'FaceAlpha',0.5);
%}

%% Plot probability distribution

[X,Y,Z] = meshgrid(-L/2:L/20:L/2, -L/2:L/20:L/2, -L/2:L/20:L/2);
V = interpolateSolution(results,X,Y,Z,eigenV);
V = reshape(V,size(X));
prob = abs(V).^2;

%{
%figure('Name', 'Sliced in x direction');
colormap jet;
subplot(2,2,2);
contourslice(X,Y,Z,V.^2,-0.45:0.01:0.45,[],[]);
title(sprintf('Probability Distribution, sliced in x direction for No:%s eigenvalue', num2str(eigenv)));
xlabel('x'); ylabel('y'); zlabel('z');
colorbar; view(-11,14); axis equal;

%figure('Name', 'Sliced in y direction');
colormap jet;
subplot(2,2,3);
contourslice(X,Y,Z,V.^2,[],-0.45:0.01:0.45,[]);
title(sprintf('Probability Distribution, sliced in y direction for No:%s eigenvalue', num2str(eigenv)));
xlabel('x'); ylabel('y'); zlabel('z');
colorbar; view(-11,14); axis equal;

%figure('Name', 'Sliced in z direction');
colormap jet;
subplot(2,2,4);
contourslice(X,Y,Z,V.^2,[],[],-0.45:0.01:0.45);
title(sprintf('Probability Distribution, sliced in z direction for No:%s eigenvalue', num2str(eigenv)));
xlabel('x'); ylabel('y'); zlabel('z');
colorbar; view(-11,14); axis equal;
%}

xslice = 0;        % location of y-z planes
yslice = 0;        % location of x-z planes
zslice = 0;        % location of x-y planes
figure('Name', 'Probability Distribution');
slice(X, Y, Z, prob, xslice, yslice, zslice)
xlabel('x'); ylabel('y'); zlabel('z'); colorbar; daspect([1 1 1]); axis equal;
title(sprintf('Probability Distribution for eigen energy level %s', num2str(eigenV)));