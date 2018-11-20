%................................................................

% MATLAB codes for Finite Element Analysis
% problem4.m
% antonio ferreira 2008

% clear memory
clear all

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=30e6;   A=2;    EA=E*A; 

% generation of coordinates and connectivities
numberElements=3;  
numberNodes=4;
elementNodes=[1 2;1 3;1 4];
nodeCoordinates=[ 0 0;0 120;120 120;120 0];
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
    
GDof=2*numberNodes; % GDof: total number of degrees of freedom
displacements=zeros(GDof,1);
force=zeros(GDof,1);

% applied load at node 2
force(2)=-10000.0;

% computation of the system stiffness matrix
 [stiffness]=...
    formStiffness2Dtruss(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);

% boundary conditions and solution
prescribedDof=[3:8]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% drawing displacements
us=1:2:2*numberNodes-1;
vs=2:2:2*numberNodes;
figure
L=xx(2)-xx(1);
%L=node(2,1)-node(1,1);
XX=displacements(us);YY=displacements(vs);
dispNorm=max(sqrt(XX.^2+YY.^2));
scaleFact=15000*dispNorm;
clf
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');

% stresses at elements
stresses2Dtruss(numberElements,elementNodes,...
    xx,yy,displacements,E)

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)
