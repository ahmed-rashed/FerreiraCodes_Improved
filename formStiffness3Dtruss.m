function  stiffness=formStiffness3Dtruss(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,A_vec)

stiffness=zeros(GDof); 
% computation of the system stiffness matrix
for iElement=1:numberElements
  % elementDof: element degrees of freedom (Dof)
  iNodes=elementNodes(iElement,:);
  elementDof=[3*iNodes(1)-2 3*iNodes(1)-1 3*iNodes(1)...
          3*iNodes(2)-2 3*iNodes(2)-1 3*iNodes(2)] ;
  x1=nodeCoordinates(iNodes(1),1);
  y1=nodeCoordinates(iNodes(1),2);
  z1=nodeCoordinates(iNodes(1),3);
  x2=nodeCoordinates(iNodes(2),1);
  y2=nodeCoordinates(iNodes(2),2);
  z2=nodeCoordinates(iNodes(2),3);
  L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +...
      (z2-z1)*(z2-z1));
  CXx = (x2-x1)/L;CYx = (y2-y1)/L;CZx = (z2-z1)/L;  
     
  T = [CXx*CXx CXx*CYx CXx*CZx ; CYx*CXx CYx*CYx CYx*CZx ; ...
         CZx*CXx CZx*CYx CZx*CZx]; 
  stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+E_vec(iElement)*A_vec(iElement)/L*[T -T ; -T T];
end 