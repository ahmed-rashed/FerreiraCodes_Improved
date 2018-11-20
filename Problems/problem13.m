%................................................................

% MATLAB codes for Finite Element Analysis
% problem13.m
% antonio ferreira 2008

% clear memory
clear all

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
    4 5 0;
    ];
xx=nodeCoordinates(:,1);   
yy=nodeCoordinates(:,2); 
zz=nodeCoordinates(:,3);
elementNodes=[1 5;2 6;3 7; 4 8; 5 6; 6 7; 7 8; 8 5]; 
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
    % GDof: global number of degrees of freedom
GDof=6*numberNodes; 
U=zeros(GDof,1);
force=zeros(GDof,1);
stiffness=zeros(GDof); 

%force vector
force(37)=-15;

% calculation of the system stiffness matrix
% and force vector
% stiffness matrix
[stiffness]=...
    formStiffness3Dframe(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);

% boundary conditions and solution
prescribedDof=[1:24];

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% displacements
disp('Displacements')
jj=1:GDof; format long
f=[jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)

%drawing mesh and deformed shape
U=displacements;   
clf
drawingMesh(nodeCoordinates+500*[U(1:6:6*numberNodes)...
    U(2:6:6*numberNodes) U(3:6:6*numberNodes)],...
    elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k--');
