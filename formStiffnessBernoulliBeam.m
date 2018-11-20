%................................................................

function  [stiffness,force]=...
    formStiffnessBernoulliBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,EI,P);

force=zeros(GDof,1);
stiffness=zeros(GDof); 
% calculation of the system stiffness matrix
% and force vector
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ 2*(indice(1)-1)+1 2*(indice(2)-1)...
      2*(indice(2)-1)+1 2*(indice(2)-1)+2]; 
 % ndof=length(indice);  
  % length of element
  LElem=xx(indice(2))-xx(indice(1))  ;
  ll=LElem;
  k1=EI/(LElem)^3*[12   6*LElem -12 6*LElem;
     6*LElem 4*LElem^2 -6*LElem 2*LElem^2;
     -12 -6*LElem 12 -6*LElem ;
     6*LElem 2*LElem^2 -6*LElem 4*LElem^2];

  f1=[P*LElem/2 P*LElem*LElem/12 P*LElem/2 ...
      -P*LElem*LElem/12]';
  
  % equivalent force vector
  force(elementDof)=force(elementDof)+f1;  
  
  % stiffness matrix
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+k1;
end
