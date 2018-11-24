% Timoshenko beam under buckling loads (P)
% antonio ferreira 2008
%Modified by Ahmed Rashed

clearvars

E=10e6;
nu=0.333;
L=1;
b=1;
h=0.001;
I=b*h^3/12;
kapa=5/6;
A=b*h;

P = 1; % uniform pressure
% constitutive matrix
G=E/2/(1+nu);
C=[ EI   0
    0    kapa*h*G];

% mesh
N_elements=40;  
nodeCoordinates=linspace(0,L,N_elements+1);
xx=nodeCoordinates';
x=xx.';
for i=1:size(nodeCoordinates,2)-1
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1
end
% generation of coordinates and connectivities
N_nodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*N_nodes; 

% computation of the system stiffness, Kg
[K_Assembly,K_G_Assembly]=formStiffnessBucklingTimoshenkoBeam(GDof,N_elements,elementNodes,N_nodes,xx,C,P,I,h);

% boundary conditions (CC)
fixedNodeW =find(xx==min(nodeCoordinates(:)) | xx==max(nodeCoordinates(:)));
fixedNodeTX=fixedNodeW;
prescribedDof=[fixedNodeW; fixedNodeTX+N_nodes];

% boundary conditions (SS)
% fixedNodeW =find(xx==min(nodeCoordinates(:))...
%     | xx==max(nodeCoordinates(:)));
% fixedNodeTX=[];
% prescribedDof=[fixedNodeW; fixedNodeTX+numberNodes];

% Buckling problem

activeDof=setdiff(1:GDof,prescribedDof);
[V,D]=eig(K_Assembly(activeDof,activeDof),K_G_Assembly(activeDof,activeDof));
D=diag(D);
[D,ii] = sort(D);
V = V(:,ii);

kapa=5/6;
PcrSS=pi^2*E*I/L^2/(1+pi^2*E*I/(L^2*kapa*G*A))
PcrCC=pi^2*E*I/(L/2)^2/(1+pi^2*E*I/(L^2/4*kapa*G*A))

modeNumber=4;

V1=zeros(GDof,1);
V1(activeDof,1:modeNumber)=V(:,1:modeNumber);
       
% drawing eigenmodes
drawEigenmodes1D(modeNumber,N_nodes,V1,xx,x)
