%................................................................

% MATLAB codes for Finite Element Analysis
% problem9a.m
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E=1e6;  L=10; t=L/1000;I=1*t^3/12;  EI=E*I;

% generation of coordinates and connectivities
numberElements=3;
nodeCoordinates=linspace(0,L,numberElements+1)';
xx=nodeCoordinates;L=max(nodeCoordinates);
for i=1:numberElements; 
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);

% distributed force
P=-1000;  

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
    % GDof: global number of degrees of freedom
GDof=2*numberNodes; 
%U=zeros(GDof,1);
stiffnessSpring=zeros(GDof+1);
forceSpring=zeros(GDof+1,1);

% stiffess matrix and force vector
[stiffness,force]=...
    formStiffnessBernoulliBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,EI,P);
% spring added
stiffnessSpring(1:GDof,1:GDof)=stiffness;
forceSpring(1:GDof)=force;
k=10;
stiffnessSpring([GDof-1 GDof+1],[GDof-1 GDof+1])=...
    stiffnessSpring([GDof-1 GDof+1],[GDof-1 GDof+1])+[k -k;-k k];

% boundary conditions and solution
fixedNodeU =[1]'; fixedNodeV =[2]';
prescribedDof=[fixedNodeU;fixedNodeV;GDof+1];

% solution
displacements=solution(GDof+1,prescribedDof,...
    stiffnessSpring,forceSpring);

% displacements
disp('Displacements')
jj=1:GDof+1; format
[jj' displacements]

% drawing deformed shape
U=displacements(1:2:2*numberNodes);   
plot(nodeCoordinates,U,'*')

% exact solution by Bathe (Solutions Manual of Procedures ...)
load=[L*P/3;L*P/3;L*P/6];
K=E*I/L^3*[189 -108 27;-108 135 -54;27 -54 27+k*L^3/E/I];
X=K\load


