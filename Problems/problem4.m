% MATLAB codes for Finite Element Analysis
% problem4.m
% antonio ferreira 2008
% Modified by Ahmed Rashed

clc
close all
clearvars

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E_vec=[30e6,30e6,30e6];   A_vec=[2,2,2];

% generation of coordinates and connectivities
numberElements=3;  
numberNodes=4;
elementNodes=[1 2;
              1 3;
              1 4];
nodeCoordinates=[0 0;
                0 120;
                120 120;
                120 0];

% GDof: total number of degrees of freedom
GDof=2*numberNodes;

% stiffness matrix
K=formStiffness2Dtruss(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,A_vec);

% boundary conditions
prescribedDof=3:8;

% force : force vector
force=nan(GDof,1);
force(1)=0;
force(2)=-1e4;

%displacement vector
displacements=nan(GDof,1);
displacements(prescribedDof)=0;

% solution
[displacements,force]=solution(prescribedDof,K,displacements,force)

% % drawing displacements
% us=1:2:2*numberNodes-1;
% vs=2:2:2*numberNodes;
% figure
% L=nodeCoordinates(2,1)-nodeCoordinates(1,1);
% %L=node(2,1)-node(1,1);
% XX=displacements(us);YY=displacements(vs);
% dispNorm=max(sqrt(XX.^2+YY.^2));
% scaleFact=15000*dispNorm;

% hold on
% drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes,'L2','k.-');
% drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');

% stresses at elements
sigma=stresses2Dtruss(numberElements,elementNodes,nodeCoordinates,displacements,E_vec)