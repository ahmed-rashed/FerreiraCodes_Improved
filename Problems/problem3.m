%................................................................

% MATLAB codes for Finite Element Analysis
% problem3.m
% ref: D. Logan, A first couse in the finite element method,
% third Edition, page 121, exercise P3-10
% direct stiffness method
% antonio ferreira 2008

% clear memory
clear all

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
% k: spring stiffness
E=70000;A=200;k=2000;

% generation of coordinates and connectivities
% numberElements: number of elements
numberElements=3;  
numberNodes=4;
elementNodes=[1 2; 2 3; 3 4];
nodeCoordinates=[0 2000 4000 4000];
xx=nodeCoordinates;

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
displacements=zeros(numberNodes,1);
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes); 

% applied load at node 2
force(2)=8000;

% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  L=nodeCoordinates(elementDof(2))-nodeCoordinates(elementDof(1));
  if e<3
      ea(e)=E*A/L;
  else
      ea(e)=k;
  end
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ea(e)*[1 -1;-1 1];
end 
% boundary conditions and solution
% prescribed dofs
prescribedDof=[1;4]; 
% free Dof : activeDof
activeDof=setdiff([1:numberNodes]',[prescribedDof]);

% solution
GDof=4;
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    numberNodes,prescribedDof)