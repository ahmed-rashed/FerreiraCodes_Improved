%................................................................

% MATLAB codes for Finite Element Analysis
% problem6.m
% ref: D. Logan, A first couse in the finite element method,
% third Edition, mixing trusses with springs
% antonio ferreira 2008

% clear memory
clearvars

% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=210000;   A=500;    EA=E*A;

% generation of coordinates and connectivities
nodeCoordinates=[0 0;-5000*cos(pi/4) 5000*sin(pi/4); -10000 0];
elementNodes=[ 1 2;1 3;1 4];
numberElements=size(elementNodes,1);
numberNodes=size(nodeCoordinates,1)+1; % spring added
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
GDof=2*numberNodes;
U=zeros(GDof,1);
force=zeros(GDof,1);
stiffness=zeros(GDof); 

% applied load at node 2
force(2)=-25000;

% computation of the system stiffness matrix
for e=1:numberElements-1; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));
  length_element=sqrt(xa*xa+ya*ya);
  C=xa/length_element;
  S=ya/length_element;   
    k1=EA/length_element*...
        [C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;
        -C*C -C*S C*C C*S;-C*S -S*S C*S S*S];    
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+k1;
end 
% spring stiffness in global Dof
stiffness([2 7],[2 7])= stiffness([2 7],[2 7])+2000*[1 -1;-1 1];

% boundary conditions and solution
prescribedDof=[3:8]';

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
    GDof,prescribedDof)

% stresses at elements
for e=1:numberElements-1                          
  indice=elementNodes(e,:);
  elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));
  length_element=sqrt(xa*xa+ya*ya);
  C=xa/length_element;
  S=ya/length_element;   
  sigma(e)=E/length_element* ...
      [-C  -S C S]*displacements(elementDof); 
end    
disp('stresses')
sigma'
