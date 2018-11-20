%................................................................

% MATLAB codes for Finite Element Analysis
% problem2.m
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E  = 30e6;A=1;EA=E*A; L  = 90; 

% generation of coordinates and connectivities
% numberElements: number of elements
numberElements=3;  
% generation equal spaced coordinates
nodeCoordinates=linspace(0,L,numberElements+1);
xx=nodeCoordinates;
% numberNodes: number of nodes
numberNodes=size(nodeCoordinates,2);

% elementNodes: connections at elements
ii=1:numberElements; 
elementNodes(:,1)=ii; 
elementNodes(:,2)=ii+1;

% for structure:
    % D_col: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
D_col=zeros(numberNodes,1);
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes); 
% applied load at node 2
force(2)=3000.0;

% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  nn=length(elementDof);    
  length_element=nodeCoordinates(elementDof(2))...
      -nodeCoordinates(elementDof(1));
  detJacobian=length_element/2;invJacobian=1/detJacobian;

  % central Gauss point (xi=0, weight W=2)
  [shape,naturalDerivatives]=shapeFunctionL2(0.0);
     Xderivatives=naturalDerivatives*invJacobian;

% B matrix     
B=zeros(1,nn);  B(1:nn)  = Xderivatives(:); 
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+B'*B*2*detJacobian*EA;
end 

% boundary conditions and solution
% prescribed dofs
fixedDof=find(xx==min(nodeCoordinates(:)) ...
    | xx==max(nodeCoordinates(:)))'; 
prescribedDof=[fixedDof]
% free Dof : activeDof
activeDof=setdiff([1:numberNodes]',[prescribedDof]);

% solution
GDof=numberNodes;
D_col=solution(GDof,prescribedDof,stiffness,force);

% output D_col/reactions
outputDisplacementsReactions(D_col,stiffness,...
    numberNodes,prescribedDof)
