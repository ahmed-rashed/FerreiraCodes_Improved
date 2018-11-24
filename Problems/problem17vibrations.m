%................................................................

% MATLAB codes for Finite Element Analysis
% problem17.m
% 2D proble: thin plate in tension
% antonio ferreira 2008

% clear memory
clearvars;colordef white;clf

% materials
L  = 1;     % comprimento 
rho=1;
E  = 10e7;     poisson = 0.30;rho=1;thickness=0.1;  
I=thickness^3/12;% momento de inercia
A=1*thickness;

% matriz C
C=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% load
P = 1e6;

%Mesh generation
Lx=L;
Ly=thickness;
numberElementsX=20;
numberElementsY=4;
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(Lx,Ly,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% calculation of the system stiffness matrix
[K_Assembly,M_Assembly]=formStiffness2D(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,C,rho,thickness);

% boundary conditions 
fixedNodeX=find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY=find(nodeCoordinates(:,1)==0);  % fixed in YY
prescribedDof=[fixedNodeX; fixedNodeY+numberNodes];
activeDof=setdiff(1:GDof,prescribedDof);

% perform eigenproblem
[V1,D1] = eig(K_Assembly(activeDof,activeDof),M_Assembly(activeDof,activeDof)); 
D1 = diag(sqrt(D1)*L*L*sqrt(rho*A/E/I));
% drawing eigenmodes
numberOfModes=4;
% sort out eigenvalues
[D1,ii] = sort(D1); ii = ii(1:numberOfModes); VV = V1(:,ii);