function  [K_asembly_local,F_equiv_local]=formStiffnessBernoulliBeam(GDof,elementNodes,x_col,E,I,P)

N_elements=size(elementNodes,1);

K_asembly_local=zeros(GDof);
F_equiv_local=zeros(GDof,1);
for iElement=1:N_elements
  i_nodes=elementNodes(iElement,:);
  elementDof=[2*(i_nodes(1)-1)+1 2*(i_nodes(2)-1) 2*(i_nodes(2)-1)+1 2*(i_nodes(2)-1)+2];

  L=x_col(i_nodes(2))-x_col(i_nodes(1));
  
  k_local=E*I/L^3*[ 12   6*L -12 6*L;
                    6*L 4*L^2 -6*L 2*L^2;
                    -12 -6*L 12 -6*L ;
                    6*L 2*L^2 -6*L 4*L^2];

  % Assembly stiffness matrix
  K_asembly_local(elementDof,elementDof)=K_asembly_local(elementDof,elementDof)+k_local;

  f_eq_1ocal=[ P*L/2
            P*L*L/12
            P*L/2
            -P*L*L/12];
  
  % equivalent force vector
  F_equiv_local(elementDof)=F_equiv_local(elementDof)+f_eq_1ocal;  
end
