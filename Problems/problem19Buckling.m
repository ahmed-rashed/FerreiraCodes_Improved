% Buckling analysis of Mindlin plates
% antonio ferreira 2008
% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

clc
clearvars
close all

E=10920;
nu=0.30;

% L: side lenght
L=1;     

h=0.001;
I=h^3/12;

% kapa: shear correction factor
kapa=5/6;         

% constitutive matrix
% bending part
D_b=I*E/(1-nu^2)*...
    [   1                   nu          0 ; 
        nu        1                     0 ; 
         0                  0    (1-nu)/2 ];

% shear part
G=E/2/(1+nu);
D_s=kapa*h*G*eye(2);

% initial stress matrix
sigmaX=1/h;
sigmaXY=0;
sigmaY=0;
sigmaMatrix=[   sigmaX sigmaXY
                sigmaXY sigmaY];
                         
% mesh generation ...
numberElementsX=10;
numberElementsY=10;
% number of elements
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes]=rectangularMesh(L, L, numberElementsX, numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
figure
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off

numberNodes=size(xx,1);    % number of nodes
GDof=3*numberNodes;        % total number of DOFs

% stiffness and geometric stiffness matrices
K_Assembly=formMatricesMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,D_s,D_b,h,I);
K_G_Assembly=formGeometricStiffnessMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,h);

% Essential boundary conditions    
[prescribedDof,activeDof,fixedNodeW]=EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);
  
% buckling analysis ...

% perform eigenproblem
[V1,D1]=eig(K_Assembly(activeDof,activeDof),...
    K_G_Assembly(activeDof,activeDof)); 
D1=diag(D1);
% drawing eigenmodes
numberOfModes=12;
% sort out eigenvalues
[D1,ii]=sort(D1); ii=ii(1:numberOfModes); 
VV=V1(:,ii);
activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
NNN=size(activeDofW);

% normalize results
disp('D1(1)/pi/pi/C_bending(1,1)')
D1(1)/pi/pi/D_b(1,1)
D1(1)*pi*pi*D_b(1,1)

VVV(1:numberNodes,1:numberOfModes)=0;
for i=1:numberOfModes
    VVV(activeDofW,i)=VV(1:NNN,i);
end
%   
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);

D1=D1/pi/pi/D_b(1,1);
% drawing Eigenmodes
drawEigenmodes2D(x,y,VVV,NN,N,D1)
