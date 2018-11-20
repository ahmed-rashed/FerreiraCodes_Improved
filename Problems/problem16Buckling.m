%................................................................

% MATLAB codes for Finite Element Analysis
% problem16Buckling.m
% Timoshenko beam under buckling loads (P)
% antonio ferreira 2008

% clear memory
clear all

% E; modulus of elasticity
% G; shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
E=10e6; poisson = 0.333;L  = 1;thickness=0.001;
I=thickness^3/12;
EI=E*I;
kapa=5/6;
A=1*thickness;
% 

P = 1; % uniform pressure
% constitutive matrix
G=E/2/(1+poisson);
C=[   EI   0; 0    kapa*thickness*G];

% mesh
numberElements     = 40;  
nodeCoordinates=linspace(0,L,numberElements+1);
xx=nodeCoordinates';x=xx';
for i=1:size(nodeCoordinates,2)-1
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1
end
% generation of coordinates and connectivities
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% computation of the system stiffness, Kg
[stiffness,Kg]=...
    formStiffnessBucklingTimoshenkoBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,C,P,I,thickness);

% boundary conditions (CC)
fixedNodeW =find(xx==min(nodeCoordinates(:))...
    | xx==max(nodeCoordinates(:)));
fixedNodeTX=fixedNodeW;
prescribedDof=[fixedNodeW; fixedNodeTX+numberNodes];

% boundary conditions (SS)
% fixedNodeW =find(xx==min(nodeCoordinates(:))...
%     | xx==max(nodeCoordinates(:)));
% fixedNodeTX=[];
% prescribedDof=[fixedNodeW; fixedNodeTX+numberNodes];

% Buckling problem

    activeDof=setdiff([1:GDof]',[prescribedDof]);
    [V,D]=eig(stiffness(activeDof,activeDof),Kg(activeDof,activeDof));
    D=diag(D);[D,ii] = sort(D);  V = V(:,ii);

kapa=5/6;
PcrSS=pi*pi*E*I/L^2*(1/(1+pi*pi*E*I/(L*L*kapa*G*A)))
PcrCC=pi*pi*E*I/(L/2)^2*(1/(1+pi*pi*E*I/(L*L/4*kapa*G*A)))
    
    modeNumber=4;

    V1=zeros(GDof,1);
    V1(activeDof,1:modeNumber)=V(:,1:modeNumber);
       
% drawing eigenmodes
drawEigenmodes1D(modeNumber,numberNodes,V1,xx,x)