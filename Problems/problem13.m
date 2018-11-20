%................................................................

% MATLAB codes for Finite Element Analysis
% problem13.m
% antonio ferreira 2008

% clear memory
clc
clearvars

% E; modulus of elasticity
% I: second moments of area
% J: polar moment of inertia
% G: shear modulus
% L: length of bar
E=210e6; A=0.02;  
Iy=10e-5;   Iz=20e-5; J=5e-5; G=84e6;

% generation of coordinates and connectivities
nodeCoordinates=[0 0 0; 
                0 0 4; 
                4 0 4; 
                4 0 0;
                0 5 0; 
                0 5 4; 
                4 5 4; 
                4 5 0];
elementNodes=[1 5;2 6;3 7; 4 8; 5 6; 6 7; 7 8; 8 5]; 
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);

% GDof: global number of degrees of freedom
GDof=6*numberNodes; 

% calculation of the system stiffness matrix
stiffness=formStiffness3Dframe(GDof,numberElements,elementNodes,nodeCoordinates,E,A,Iz,Iy,G,J);

% boundary conditions and solution
prescribedDof=1:4*6;

%force vector
force=nan(GDof,1);
force(4*6+1:end)=0;
force(37)=-15;

%displacement vector
displacements=nan(GDof,1);
displacements(prescribedDof)=0;

% solution
[displacements,force]=solution(GDof,prescribedDof,stiffness,displacements,force)

%drawing mesh and deformed shape
clf
drawingMesh(nodeCoordinates+500*[displacements(1:6:6*numberNodes)...
    displacements(2:6:6*numberNodes) displacements(3:6:6*numberNodes)],...
    elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k--');
