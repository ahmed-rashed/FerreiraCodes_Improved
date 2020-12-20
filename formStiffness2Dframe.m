function  K_assembly=formStiffness2Dframe(GDof,numberElements,elementNodes,nodeCoordinates,E_vec,I_vec,A_vec)
%This corrects the strange node ordering of Ferreira's book

N_node_DOF=3;   %2D problem
K_assembly=zeros(GDof); 
for iElement=1:numberElements 
    iNodes=elementNodes(iElement,:);
    x_local=nodeCoordinates(iNodes(2),:)-nodeCoordinates(iNodes(1),:);
    L=norm(x_local);
  
    E=E_vec(iElement);
    A=A_vec(iElement);
    II=I_vec(iElement);
    k_local=E*[   A/L,0,0,-A/L,0,0
                0,12*II/L^3,6*II/L^2,0,-12*II/L^3,6*II/L^2
                0,6*II/L^2,4*II/L,0,-6*II/L^2,2*II/L
                -A/L,0,0,A/L,0,0
                0,-12*II/L^3,-6*II/L^2,0,12*II/L^3,-6*II/L^2
                0,6*II/L^2,2*II/L,0,-6*II/L^2,4*II/L];
    l=x_local(1)/L;
    m=x_local(2)/L;

    T_node=[l,m,0
            -m,l,0
            0,0,1];
    T_element=[T_node,zeros(N_node_DOF);zeros(N_node_DOF),T_node];

    k_global=T_element.'*k_local*T_element;
    
    elementDof=[N_node_DOF*iNodes(1)+(-N_node_DOF:-1)+1,N_node_DOF*iNodes(2)+(-N_node_DOF:-1)+1];
    K_assembly(elementDof,elementDof)=K_assembly(elementDof,elementDof)+k_global;
end
