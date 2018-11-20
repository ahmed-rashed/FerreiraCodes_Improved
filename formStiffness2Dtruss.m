function K=formStiffness2Dtruss(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,A_vec)
% Computes the assembly stiffness matrix

K=zeros(GDof,GDof);
for iElement=1:numberElements
  iNodes=elementNodes(iElement,:);       
  elementDof=[ iNodes(1)*2-1 iNodes(1)*2 iNodes(2)*2-1 iNodes(2)*2] ;

  D_x=nodeCoordinates(iNodes(2),1)-nodeCoordinates(iNodes(1),1);
  D_y=nodeCoordinates(iNodes(2),2)-nodeCoordinates(iNodes(1),2);
  L=sqrt(D_x^2+D_y^2);
  l=D_x/L;
  m=D_y/L;
  k=E_vec(iElement)*A_vec(iElement)/L*[l*l l*m -l*l -l*m;
            l*m m*m -l*m -m*m;
            -l*l -l*m l*l l*m;
            -l*m -m*m l*m m*m];
  K(elementDof,elementDof)=K(elementDof,elementDof)+k;
end 