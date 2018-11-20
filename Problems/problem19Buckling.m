%................................................................

% MATLAB codes for Finite Element Analysis
% problem19Buckling.m
% this function performs buckling analysis
% of Mindlin plates using 3 degrees of freedom per node
% antonio ferreira 2008

clearvars
colordef white

% material properties
% modulusOfElasticity  = Young's modulus
% PoissonRatio         = Poisson's ratio

modulusOfElasticity  = 10920;  % Young
PoissonRatio = 0.30;  % coef. Poisson

% L: side lenght
L  = 1;     

thickness=0.001;
I=thickness^3/12;

% kapa: shear correction factor
kapa=5/6;         

% constitutive matrix
% bending part
C_bending=...
    I*modulusOfElasticity/(1-PoissonRatio^2)*...
    [   1                   PoissonRatio          0 ; 
        PoissonRatio        1                     0 ; 
         0                  0    (1-PoissonRatio)/2 ];
                 
% shear part
C_shear=...
kapa*thickness*modulusOfElasticity/2/(1+PoissonRatio)*eye(2);

% initial stress matrix
sigmaX=1/thickness;
sigmaXY=0;
sigmaY=0;
sigmaMatrix=[ sigmaX sigmaXY; sigmaXY sigmaY];
                         
% mesh generation ...
% numberElementsX: number of elements in x
% numberElementsY: number of elements in y
numberElementsX=10;
numberElementsY=10;
% number of elements
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L, L, numberElementsX, numberElementsY);
xx=nodeCoordinates(:,1);   yy=nodeCoordinates(:,2);
figure
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off

numberNodes=size(xx,1);    % number of nodes
GDof=3*numberNodes;        % total number of DOFs

% stiffness and geometric stiffness matrices
[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I);

[geometric]=...
formGeometricStiffnessMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness);

% Essential boundary conditions    
[prescribedDof,activeDof,fixedNodeW]=...
EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);
  
% buckling analysis ...

% perform eigenproblem
[V1,D1] = eig(stiffness(activeDof,activeDof),...
    geometric(activeDof,activeDof)); 
D1 = diag(D1);
% drawing eigenmodes
numberOfModes=12;
% sort out eigenvalues
[D1,ii] = sort(D1); ii = ii(1:numberOfModes); 
VV = V1(:,ii);
activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
NNN=size(activeDofW);

% normalize results
disp('D1(1)/pi/pi/C_bending(1,1)')
D1(1)/pi/pi/C_bending(1,1)
D1(1)*pi*pi*C_bending(1,1)

VVV(1:numberNodes,1:numberOfModes)=0;
for i=1:numberOfModes
    VVV(activeDofW,i)=VV(1:NNN,i);
end
%   
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);

D1=D1/pi/pi/C_bending(1,1);
% drawing Eigenmodes
drawEigenmodes2D(x,y,VVV,NN,N,D1)
