%................................................................

% MATLAB codes for Finite Element Analysis
% problem8.m
% ref: D. Logan, A first couse in the finite element method,
% third Edition, A second 3D truss example
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=210000; 
A=[100 100 100 100]; % area for various sections

% generation of coordinates and connectivities
nodeCoordinates=[4000 4000 3000;
                    0 4000 0;
                    0 4000 6000;
                 4000 0    3000;
                 8000 -1000 1000];
elementNodes=[1 2;1 3;1 4;1 5]; 
numberElements=size(elementNodes,1);
numberNodes=size(nodeCoordinates,1); 
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
    % GDof: global number of degrees of freedom
GDof=3*numberNodes; 
U=zeros(GDof,1);
force=zeros(GDof,1);

% applied load at node 2
force(2)=-10000;

% stiffness matrix
[stiffness]=...
    formStiffness3Dtruss(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A);

% boundary conditions and solution
prescribedDof=[4:15]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% stresses at elements
stresses3Dtruss(numberElements,elementNodes,nodeCoordinates,...
    displacements,E)
