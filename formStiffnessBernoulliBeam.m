function  [K_asembly_local,F_equivalent_local]=formStiffnessBernoulliBeam(GDof,N_elements,elementNodes,xx,E,I,P)

K_asembly_local=zeros(GDof);
F_equivalent_local=zeros(GDof,1);
for i_element=1:N_elements
  i_nodes=elementNodes(i_element,:);
  elementDof=[2*(i_nodes(1)-1)+1 2*(i_nodes(2)-1) 2*(i_nodes(2)-1)+1 2*(i_nodes(2)-1)+2];

  L=xx(i_nodes(2))-xx(i_nodes(1));
  
  k_local=E*I/L^3*[ 12   6*L -12 6*L;
                    6*L 4*L^2 -6*L 2*L^2;
                    -12 -6*L 12 -6*L ;
                    6*L 2*L^2 -6*L 4*L^2];

  % Assembly stiffness matrix
  K_asembly_local(elementDof,elementDof)=K_asembly_local(elementDof,elementDof)+k_local;

  f_equivalent_1ocal=[ P*L/2
            P*L*L/12
            P*L/2
            -P*L*L/12];
  
  % equivalent force vector
  F_equivalent_local(elementDof)=F_equivalent_local(elementDof)+f_equivalent_1ocal;  
end
