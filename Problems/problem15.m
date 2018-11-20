%................................................................

% MATLAB codes for Finite Element Analysis
% problem15.m
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% I: second moments of area
% J: polar moment of inertia
% G: shear modulus
% L: length of bar
E=210e6; G=84e6;  I=20e-5;   J=5e-5;

% generation of coordinates and connectivities
nodeCoordinates=[4 4; 0 4; 0 0 ; 4 0];
xx=nodeCoordinates(:,1);   
yy=nodeCoordinates(:,2); 
elementNodes=[1 2; 3 1; 4 1]; 
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 

force=zeros(GDof,1);
%force vector
force(1)=-20;

% computation of the system stiffness matrix
stiffness=formStiffnessGrid(GDof,numberElements,...
    elementNodes,xx,yy,E,I,G,J);

% boundary conditions
prescribedDof=[4:12]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force)

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% % forces in elements
disp('forces in elements ')
EF=forcesInElementGrid(numberElements,elementNodes,...
    xx,yy,E,I,G,J,displacements)
