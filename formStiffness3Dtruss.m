function  [stiffness]=...
    formStiffness3Dtruss(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A);

stiffness=zeros(GDof); 
% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[3*indice(1)-2 3*indice(1)-1 3*indice(1)...
          3*indice(2)-2 3*indice(2)-1 3*indice(2)] ;
  x1=nodeCoordinates(indice(1),1);
  y1=nodeCoordinates(indice(1),2);
  z1=nodeCoordinates(indice(1),3);
  x2=nodeCoordinates(indice(2),1);
  y2=nodeCoordinates(indice(2),2);
  z2=nodeCoordinates(indice(2),3);
  L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +...
      (z2-z1)*(z2-z1));
  CXx = (x2-x1)/L;CYx = (y2-y1)/L;CZx = (z2-z1)/L;  
     
  T = [CXx*CXx CXx*CYx CXx*CZx ; CYx*CXx CYx*CYx CYx*CZx ; ...
         CZx*CXx CZx*CYx CZx*CZx]; 
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+E*A(e)/L*[T -T ; -T T];
end 
