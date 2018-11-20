function  stiffness=formStiffness2Dframe(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,I_vec,A_vec)

N_node_DOF=3;   %2D problem
stiffness=zeros(GDof); 
for iElement=1:numberElements 
    iNodes=elementNodes(iElement,:);
    D_x=nodeCoordinates(iNodes(2),1)-nodeCoordinates(iNodes(1),1);
    D_y=nodeCoordinates(iNodes(2),2)-nodeCoordinates(iNodes(1),2);
    L=sqrt(D_x^2+D_y^2);

    oneu=[1 -1;-1 1];
    oneu2=[1 -1;1 -1];
    oneu3=[1 1;-1 -1];
    oneu4=[4 2;2 4];

    k_local=[E_vec(iElement)*A_vec(iElement)/L*oneu   zeros(2,4);
    zeros(2) 12*E_vec(iElement)*I_vec(iElement)/L^3*oneu 6*E_vec(iElement)*I_vec(iElement)/L^2*oneu3;
    zeros(2) 6*E_vec(iElement)*I_vec(iElement)/L^2*oneu2 E_vec(iElement)*I_vec(iElement)/L*oneu4];

    l=D_x/L;
    m=D_y/L;

    T= [l*eye(2) m*eye(2) zeros(2);
    -m*eye(2) l*eye(2) zeros(2);
    zeros(2,4) eye(2)];

    k_global=T.'*k_local*T;
    
    elementDof=[N_node_DOF*iNodes(1)+(-N_node_DOF:-1)+1,N_node_DOF*iNodes(2)+(-N_node_DOF:-1)+1];
    stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+k_global;
end
