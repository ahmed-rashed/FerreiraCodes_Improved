%................................................................

% MATLAB codes for Finite Element Analysis
% problem1.m
% antonio ferreira 2008

% clear memory
clear all

% elementNodes: connections at elements
elementNodes=[1 2;2 3;2 4];

% numberElements: number of Elements
numberElements=size(elementNodes,1); 

% numberNodes: number of nodes
numberNodes=4;

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
displacements=zeros(numberNodes,1);
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes); 

% applied load at node 2
force(2)=10.0;

% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+[1 -1;-1 1];
end 

% boundary conditions and solution
% prescribed dofs
prescribedDof=[1;3;4]; 
% free Dof : activeDof
activeDof=setdiff([1:numberNodes]',[prescribedDof]);

% solution
displacements=stiffness(activeDof,activeDof)\force(activeDof);

% positioning all displacements
displacements1=zeros(numberNodes,1);
displacements1(activeDof)=displacements;

% output displacements/reactions
outputDisplacementsReactions(displacements1,stiffness,...
    numberNodes,prescribedDof)