%................................................................

% MATLAB codes for Finite Element Analysis
% problem19.m
% Mindlin plate in bending
% antonio ferreira 2008

% clear memory
clearvars;colordef white;clf

% materials
E  = 10920;     poisson = 0.30; kapa=5/6; 
thickness=0.1;
I=thickness^3/12;

% matrix C
% bending part
C_bending=I*E/(1-poisson^2)*...
    [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];           
% shear part
C_shear=kapa*thickness*E/2/(1+poisson)*eye(2);

% load
P = -1;

%Mesh generation
L  = 1;    
numberElementsX=20;
numberElementsY=20;
numberElements=numberElementsX*numberElementsY;
%
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 

% computation of the system stiffness matrix and force vector
[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I);

[force]=...
    formForceVectorMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P);

% % boundary conditions 
[prescribedDof,activeDof]=...
    EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
D_col=solution(GDof,prescribedDof,stiffness,force);

% D_col
disp('Displacements')
jj=1:GDof; format
f=[jj; D_col'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)

% deformed shape
figure
plot3(xx,yy,D_col(1:numberNodes),'.')
format long
D1=E*thickness^3/12/(1-poisson^2);
min(D_col(1:numberNodes))*D1/L^4
