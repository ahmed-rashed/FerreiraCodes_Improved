%................................................................

% MATLAB codes for Finite Element Analysis
% problem1.m
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
    % D_col: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
D_col=zeros(numberNodes,1);
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
D_col=stiffness(activeDof,activeDof)\force(activeDof);

% positioning all D_col
D_col1=zeros(numberNodes,1);
D_col1(activeDof)=D_col;

% output D_col/reactions
outputDisplacementsReactions(D_col1,stiffness,...
    numberNodes,prescribedDof)
