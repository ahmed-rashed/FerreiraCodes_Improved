% Mindlin plate in free vibrations
% antonio ferreira 2008
% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

clc
clearvars
close all

% materials
E=10920;
nu=0.30;  
h=0.1;
I=h^3/12;
rho=1;
% kapa=0.8601; % cccc / cccf case
% kapa=0.822; % scsc case
kapa=5/6;  % ssss case 

% matrix D
% bending part
D_b=I*E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];           

% shear part
G=E/2/(1+nu);
D_s=kapa*h*G*eye(2);

%Mesh generation
L=1;    
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes]=rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 

% computation of the system stiffness and mass matrices
K_Assembly=formMatricesMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,D_s,D_b,h,I);

M_Assembly=formMassMatrixMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,h,rho,I);

% % boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);

G=E/2.6;
% V : mode shape
% D : frequency
% 
numberOfModes=12;
[V,D]=eig(K_Assembly(activeDof,activeDof),M_Assembly(activeDof,activeDof)); 
D=diag(sqrt(D)*L*sqrt(rho/G));
[D,ii]=sort(D); ii=ii(1:numberOfModes); 
VV=V(:,ii);
activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
NNN=size(activeDofW);
    
    VVV(1:numberNodes,1:12)=0;
    for i=1:12
        VVV(activeDofW,i)=VV(1:NNN,i);
    end
    
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);

% drawing Eigenmodes
drawEigenmodes2D(x,y,VVV,NN,N,D)
