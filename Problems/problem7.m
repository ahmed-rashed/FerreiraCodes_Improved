%................................................................

% MATLAB codes for Finite Element Analysis
% problem7.m
% ref: D. Logan, A first couse in the finite element method,
% third Edition, A 3D truss example
% antonio ferreira 2008

% clear memory
clear all

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=1.2e6; 
A=[0.302;0.729;0.187]; % area for various sections

% generation of coordinates and connectivities
nodeCoordinates=[72 0 0; 0 36 0;  0 36 72; 0 0 -48];
elementNodes=[1 2;1 3;1 4]; 
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
force(3)=-1000;

% stiffness matrix
[stiffness]=...
    formStiffness3Dtruss(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A);

% boundary conditions and solution
prescribedDof=[2 4:12]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% stresses at elements
stresses3Dtruss(numberElements,elementNodes,nodeCoordinates,...
    displacements,E)