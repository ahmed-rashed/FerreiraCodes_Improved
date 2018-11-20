% MATLAB codes for Finite Element Analysis
% antonio ferreira 2008
%Modified by Ahmed Rashed

clc
clearvars

E_vec=210e6*[1,1,1];
A_vec=0.02*[1,1,1];
Iy_vec=10e-5*[1,1,1];
Iz_vec=20e-5*[1,1,1];
J_vec=5e-5*[1,1,1];
G_vec=84e6*[1,1,1];

% generation of coordinates and connectivities
nodeCoordinates=[0 0 0;
    3 0 0;
    0 0 -3;
    0 -4 0];
elementNodes=[1 2;
                1 3;
                1 4];
numberNodes=size(nodeCoordinates,1);

% GDof: global number of degrees of freedom
GDof=6*numberNodes;

% Assembly stiffness matrix
K=formStiffness3Dframe(GDof,size(elementNodes,1),elementNodes,nodeCoordinates,E_vec,A_vec,Iz_vec,Iy_vec,G_vec,J_vec);

% boundary conditions
prescribedDof=7:24;

%force vector
force=nan(GDof,1);
force([2,4:6])=0;
force(1)=-10;
force(3)=20;

%displacement vector
displacements=nan(GDof,1);
displacements(prescribedDof)=0;

% solution
[displacements,force]=solution(prescribedDof,K,displacements,force)