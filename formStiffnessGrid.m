%................................................................

function stiffness=formStiffnessGrid(GDof,...
    numberElements,elementNodes,xx,yy,E,I,G,J)

% function to form global stiffness for grid element
stiffness=zeros(GDof); 
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[...
      (indice(1)-1)*3+1 (indice(1)-1)*3+2 (indice(1)-1)*3+3 ...
      (indice(2)-1)*3+1 (indice(2)-1)*3+2 (indice(2)-1)*3+3] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));  
  L=sqrt(xa*xa+ya*ya);
  C=xa/L;
  S=ya/L;

a1 = 12*E*I/(L*L*L);    
a2 = 6*E*I/(L*L);
a3 = G*J/L; 
a4 = 4*E*I/L;
a5 = 2*E*I/L;
% stiffness in local axes
k = [a1 0 a2 -a1 0 a2 ; 0 a3 0 0 -a3 0 ;
     a2 0 a4 -a2 0 a5 ; -a1 0 -a2 a1 0 -a2 ;
     0 -a3 0 0 a3 0; a2 0 a5 -a2 0 a4];

% transformation matrix
a=[1 0 0; 0 C S;0 -S C];
R=[a zeros(3);zeros(3) a];

  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+R'*k*R;
end 