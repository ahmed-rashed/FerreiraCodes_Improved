% Timoshenko beam in bending
% antonio ferreira 2008
%Modified by Ahmed Rashed

clearvars

E=10e7;
poisson = 0.30;
L  = 1;
thickness=0.001;
I=thickness^3/12;
kapa=5/6;
 

P=-1; % uniform distribute load

% constitutive matrix
G=E/2/(1+poisson);
C=[E*I   0;
    0    kapa*thickness*G];

% mesh
numberNodes=101;
xx=linspace(0,L,numberNodes).';

numberElements=numberNodes-1;
elementNodes=nan(numberElements,2);
elementNodes(:,1)=1:numberElements;
elementNodes(:,2)=2:numberElements+1;

% generation of coordinates and connectivities


% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% computation of the system stiffness matrix 
[stiffness,force]=...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,C,P,1,I,thickness);

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
D_col=solution(GDof,prescribedDof,stiffness,force);

% output D_col/reactions
outputDisplacementsReactions(D_col,stiffness,...
    GDof,prescribedDof)

U=D_col;
ws=1:numberNodes;

% max displacement
disp(' max displacement')
min(U(ws))

