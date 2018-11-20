%................................................................

% MATLAB codes for Finite Element Analysis
% problem16vibrations.m
% Timoshenko beam in free vibrations
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% G; shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
E=10e7; poisson = 0.30;L  = 1;thickness=0.001;
I=thickness^3/12;
EI=E*I;
kapa=5/6;
rho=1;
A=1*thickness;
% 

P = -1; % uniform pressure
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

% computation of the system stiffness, force, mass
[stiffness,force,mass]=...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,C,P,rho,I,thickness);

% boundary conditions (simply-supported at both bords)
%fixedNodeW =[1 ; numberNodes];
%fixedNodeTX=[]; 
% boundary conditions (clamped at both bords)
fixedNodeW =[1 ; numberNodes];
fixedNodeTX=fixedNodeW; 
% boundary conditions (cantilever)
fixedNodeW =[1];
fixedNodeTX=[1];; 
prescribedDof=[fixedNodeW; fixedNodeTX+numberNodes];

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% free vibration problem
activeDof=setdiff([1:GDof]',[prescribedDof]);
modeNumber=4;

[V,D]=eig(stiffness(activeDof,activeDof),...
    mass(activeDof,activeDof));

D = diag(sqrt(D)*L*L*sqrt(rho*A/E/I));
[D,ii] = sort(D); 

V1=zeros(GDof,1);
V1(activeDof,1:modeNumber)=V(:,1:modeNumber);

% drawing eigenmodes
drawEigenmodes1D(modeNumber,numberNodes,V1,xx,x)


    

   



