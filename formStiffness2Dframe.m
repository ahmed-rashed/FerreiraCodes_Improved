function  [stiffness]=...
    formStiffness2Dframe(GDof,numberElements,...
    elementNodes,numberNodes,xx,yy,EI,EA);

stiffness=zeros(GDof); 
% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ indice indice+numberNodes indice+2*numberNodes] ;
  nn=length(indice);  
  xa=xx(indice(2))-xx(indice(1))
  ya=yy(indice(2))-yy(indice(1));  
  length_element=sqrt(xa*xa+ya*ya);
  cosa=xa/length_element;  
  sena=ya/length_element;
  ll=length_element;
  
 L= [cosa*eye(2) sena*eye(2) zeros(2);
     -sena*eye(2) cosa*eye(2) zeros(2);
     zeros(2,4) eye(2)];
          
oneu=[1 -1;-1 1];
oneu2=[1 -1;1 -1];
oneu3=[1 1;-1 -1];
oneu4=[4 2;2 4];

k1=[EA/ll*oneu   zeros(2,4);
    zeros(2) 12*EI/ll^3*oneu 6*EI/ll^2*oneu3;
    zeros(2) 6*EI/ll^2*oneu2 EI/ll*oneu4];

    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+L'*k1*L;
end 

