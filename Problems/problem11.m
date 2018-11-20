%................................................................

% MATLAB codes for Finite Element Analysis
% problem11.m
% 2D frame
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E=210000; A=200; I=2e8; EA=E*A; EI=E*I;

% generation of coordinates and connectivities
numberElements=3;
nodeCoordinates=[0 0;0 6000;6000 6000;6000 0];
xx=nodeCoordinates;
for i=1:numberElements; 
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end
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

%force vector
force(2)=15000;
force(10)=10e6;

% stiffness matrix
[stiffness]=...
    formStiffness2Dframe(GDof,numberElements,...
    elementNodes,numberNodes,xx,yy,EI,EA);

% boundary conditions and solution
prescribedDof=[1 4 5 8 9 12]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

%drawing mesh and deformed shape
clf
U=displacements;
drawingMesh(nodeCoordinates+500*[U(1:numberNodes)...
    U(numberNodes+1:2*numberNodes)],elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k--');


