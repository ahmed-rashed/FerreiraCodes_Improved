% Mindlin plate in bending
% antonio ferreira 2008
% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

clc
clearvars
close all

% materials
E =10920;
nu=0.30;
kappa=5/6; 

h=0.1;

% matrix D
% bending part
D_b=E/(1-nu^2)*[1 nu 0
                nu 1 0
                0 0 (1-nu)/2];           
% shear part
G=E/2/(1+nu);
D_s=kappa*h*G*eye(2);

% load
P=-1;

%Mesh generation
L=1;    
N_elements_X=20;
N_elements_Y=20;
N_elements=N_elements_X*N_elements_Y;

[nodeCoordinates,elementNodes]=rectangularMesh(L,L,N_elements_X,N_elements_Y);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
N_nodes=size(nodeCoordinates,1);

% GDof: global number of degrees of freedom
GDof=3*N_nodes; 

[K_Assembly,M_Assembly,F_equiv]=formMatricesMindlinQ4(GDof,elementNodes,nodeCoordinates,h,D_s,D_b,P);

% boundary conditions
[prescribedDof,activeDof]=EssentialBC('ssss',GDof,nodeCoordinates(:,1),nodeCoordinates(:,2),nodeCoordinates,N_nodes);

% solution
D_col=solution(prescribedDof,K_Assembly,F_equiv);

% D_col
disp('Displacements')
jj=1:GDof;
f=[jj; D_col'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)

% deformed shape
figure
plot3(nodeCoordinates(:,1),nodeCoordinates(:,2),D_col(1:N_nodes),'.')
D1=E*h^3/12/(1-nu^2);
min(D_col(1:N_nodes))*D1/L^4
