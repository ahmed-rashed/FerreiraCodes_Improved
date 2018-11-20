%................................................................

% MATLAB codes for Finite Element Analysis
% problem12.m
% antonio ferreira 2008

% clear memory
clear all

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E=210e6; A=0.02;  
Iy=10e-5;   Iz=20e-5; J=5e-5; G=84e6;

% generation of coordinates and connectivities
nodeCoordinates=[0 0 0;    3 0 0;  0 0 -3; 0 -4 0];
xx=nodeCoordinates(:,1);   
yy=nodeCoordinates(:,2); 
zz=nodeCoordinates(:,3);
elementNodes=[1 2;1 3;1 4]; 
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
stiffness=zeros(GDof,GDof); 

%force vector
force(1)=-10;
force(3)=20;

% stiffness matrix
[stiffness]=...
    formStiffness3Dframe(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);

% boundary conditions and solution
prescribedDof=[7:24];

% solution
displacements=solution(GDof,prescribedDof,stiffness,force)

% displacements
disp('Displacements')
jj=1:GDof; format long
f=[jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)