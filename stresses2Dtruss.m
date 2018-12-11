function stress=stresses2Dtruss(N_elements,elementNodes,nodesCoordinates,D_col,E_vec)

stress=nan(N_elements,1);
for iElement=1:N_elements
  iNodes=elementNodes(iElement,:);
  elementDofs=[iNodes(1)*2-1 iNodes(1)*2 iNodes(2)*2-1 iNodes(2)*2] ;
  xa=nodesCoordinates(iNodes(2),1)-nodesCoordinates(iNodes(1),1);
  ya=nodesCoordinates(iNodes(2),2)-nodesCoordinates(iNodes(1),2);
  L=sqrt(xa*xa+ya*ya);
  C=xa/L;
  S=ya/L;   
  stress(iElement)=E_vec(iElement)/L*[-C  -S C S]*D_col(elementDofs); 
end
