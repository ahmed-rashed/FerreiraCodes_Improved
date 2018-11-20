function sigma=stresses2Dtruss(numberElements,elementNodes,nodeCoordinates,displacements,E_vec)

sigma=nan(numberElements,1);
for e=1:numberElements
  iNodes=elementNodes(e,:);
  elementDof=[iNodes(1)*2-1 iNodes(1)*2 iNodes(2)*2-1 iNodes(2)*2] ;
  xa=nodeCoordinates(iNodes(2),1)-nodeCoordinates(iNodes(1),1);
  ya=nodeCoordinates(iNodes(2),2)-nodeCoordinates(iNodes(1),2);
  L=sqrt(xa*xa+ya*ya);
  C=xa/L;
  S=ya/L;   
  sigma(e)=E_vec(e)/L*[-C  -S C S]*displacements(elementDof); 
end