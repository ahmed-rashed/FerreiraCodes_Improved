function  K=formStiffness3Dframe(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,A_vec,Iz_vec,Iy_vec,G_vec,J_vec)
% Computes the assembly stiffness matrix

K=zeros(GDof,GDof);
for iElement=1:numberElements
    iNodes=elementNodes(iElement,:);
      
    L=norm(nodeCoordinates(iNodes(2),:)-nodeCoordinates(iNodes(1),:));
    k_local=k_local_3D_Frame(L,E_vec(iElement),A_vec(iElement),Iz_vec(iElement),Iy_vec(iElement),G_vec(iElement),J_vec(iElement));
    k_global=k_global_3D_Frame(iNodes,nodeCoordinates,k_local);

    elementDof=[6*iNodes(1)+(-5:0),6*iNodes(2)+(-5:0)];
    K(elementDof,elementDof)=K(elementDof,elementDof)+k_global;
end  

