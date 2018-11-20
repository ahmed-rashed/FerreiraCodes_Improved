% MATLAB codes for Finite Element Analysis
% antonio ferreira 2008
%Modified by Ahmed Rashed

close all
clc
clearvars

E_vec=210000*[1,1,1];
A_vec=100*[1,1,1];
I_vec=2e8*[1,1,1];

% generation of coordinates and connectivities
x1=3000*(1+cos(pi/4));
nodeCoordinates=[0 3000;
                 3000 3000;
                 x1 0;
                 x1+3000 0];
elementNodes=[1 2;
              2 3;
              3 4];
numberNodes=size(nodeCoordinates,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% Assembly stiffness matrix
K=formStiffness2Dframe(GDof,size(elementNodes,1),elementNodes,nodeCoordinates,E_vec,I_vec,A_vec);

% boundary conditions
prescribedDof=[1 5 9 4 8 12];

%force vector
force=nan(GDof,1);
force(2)=0;
force(6)=-10000;
force(10)=-5e6;
force(3)=0;
force(7)=-10000;
force(11)=5e6;

%displacement vector
displacements=nan(GDof,1);
displacements(prescribedDof)=0;

% solution
[displacements,force]=solution(prescribedDof,K,displacements,force)

% % drawing undeformed and deformed meshes
% U=displacements;  
% drawingMesh(nodeCoordinates+500*[U(1:numberNodes)...
%     U(numberNodes+1:2*numberNodes)],elementNodes,'L2','k.-');
% drawingMesh(nodeCoordinates,elementNodes,'L2','k--');