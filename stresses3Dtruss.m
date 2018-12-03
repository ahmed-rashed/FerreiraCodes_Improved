function stresses3Dtruss(N_elements,elementNodes,nodesCoordinates,D_col,E_vec)

member_stress=nan(N_elements,1);
for iElement=1:N_elements
	iNodes=elementNodes(iElement,:);
	elementDof=[3*iNodes(1)-2 3*iNodes(1)-1 3*iNodes(1) 3*iNodes(2)-2 3*iNodes(2)-1 3*iNodes(2)] ;
	x1=nodesCoordinates(iNodes(1),1);
	y1=nodesCoordinates(iNodes(1),2);
	z1=nodesCoordinates(iNodes(1),3);
	x2=nodesCoordinates(iNodes(2),1);
	y2=nodesCoordinates(iNodes(2),2);
	z2=nodesCoordinates(iNodes(2),3);
	L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +(z2-z1)*(z2-z1));
	CXx = (x2-x1)/L;
	CYx = (y2-y1)/L;
	CZx = (z2-z1)/L;
 
	u=D_col(elementDof);
	member_stress(iElement)=E_vec(iElement)/L*[-CXx -CYx -CZx CXx CYx CZx]*u;
	fprintf('%3d %12.8f\n',iElement, member_stress(iElement));
end
