%................................................................

% MATLAB codes for Finite Element Analysis
% problem1b.m
% same as problem 1, but with penalty method for
% imposition of boundary conditions
% antonio ferreira 2008

% clear memory
clearvars

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

% calculation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+[1 -1;-1 1];
end 
% stiffness1: K with penalty method
stiffness1=stiffness;

% boundary conditions and solution
% prescribed dofs
prescribedDof=[1;3;4]; 
stiffness1=stiffness;
for i=1:size(prescribedDof,1)
stiffness1(prescribedDof(i),prescribedDof(i))=...
    stiffness(prescribedDof(i),prescribedDof(i))+1e10;
end
force(prescribedDof)=0;
displacements=stiffness1\force;

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    numberNodes,prescribedDof)
