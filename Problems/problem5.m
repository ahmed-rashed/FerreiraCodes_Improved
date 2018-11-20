%................................................................

% MATLAB codes for Finite Element Analysis
% problem5.m
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=70000;   A=300;    EA=E*A;

% generation of coordinates and connectivities
elementNodes=[ 1 2;1 3;2 3;2 4;1 4;3 4;3 6;4 5;4 6;3 5;5 6];
nodeCoordinates=[ 0 0;0 3000;3000 0;3000 3000;6000 0;6000 3000];
numberElements=size(elementNodes,1);
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
GDof=2*numberNodes;
U=zeros(GDof,1);
force=zeros(GDof,1);
% applied load at node 2
force(4)=-50000;
force(8)=-100000;
force(12)=-50000;

% computation of the system stiffness matrix
 [stiffness]=...
    formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);

% boundary conditions and solution
prescribedDof=[1 2 10]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);
us=1:2:2*numberNodes-1;
vs=2:2:2*numberNodes;

% drawing displacements

figure
L=xx(2)-xx(1);
%L=node(2,1)-node(1,1);
XX=displacements(us);YY=displacements(vs);
dispNorm=max(sqrt(XX.^2+YY.^2));
scaleFact=2*dispNorm;
clf
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY],...
    elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% stresses at elements
stresses2Dtruss(numberElements,elementNodes,...
    xx,yy,displacements,E)
