%................................................................

% MATLAB codes for Finite Element Analysis
% problem20.m
% laminated plate: Srinivas problem:
% S. Srinivas, A refined analysis of composite laminates, 
% J. Sound and Vibration, 30 (1973),495--507.

% antonio ferreira 2008

% clear memory
clear all;colordef white;clf

% materials
thickness=0.1;

% load
P = -1;

%Mesh generation
L  = 1;    
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=5*numberNodes; 

% computation of the system stiffness matrix
% the shear correction factors are automatically 
% calculted for any laminate

[AMatrix,BMatrix,DMatrix,SMatrix,qbarra]=srinivasMaterial(thickness);

stiffness=formStiffnessMatrixMindlinQ45laminated5dof...
    (GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,AMatrix,...
    BMatrix,DMatrix,SMatrix);

% computation of the system force vector
[force]=...
    formForceVectorMindlinQ45dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P);

% boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=...
    EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
U=solution(GDof,prescribedDof,stiffness,force);

% drawing deformed shape and normalize results
% to compare with Srinivas
ws=1:numberNodes;
disp('maximum displacement')
abs(min(U(ws))*0.999781/thickness)
figure (1)
plot3(xx,yy,U(ws),'.')

% stress computation (Srinivas only)
disp('stress computation (Srinivas only)')
SrinivasStress(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,qbarra,U,thickness);