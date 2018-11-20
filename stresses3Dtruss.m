function stresses3Dtruss(numberElements,elementNodes,...
    nodeCoordinates,displacements,E)

% stresses in 3D truss elements
fprintf('Stresses in elements\n')
ff=zeros(numberElements,6); format
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
 
  u=displacements(elementDof);
    member_stress(e)=E/L*[-CXx -CYx -CZx CXx CYx CZx]*u;
    fprintf('%3d %12.8f\n',e, member_stress(e));
end