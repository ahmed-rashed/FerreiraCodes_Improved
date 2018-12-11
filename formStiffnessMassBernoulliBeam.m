function  [K_Assembly_local,F_equiv_local,M_Assembly_local]=formStiffnessMassBernoulliBeam(GDof,elementNodes,x_col,E,I,P,rho,A)

N_elements=size(elementNodes,1);

K_Assembly_local=zeros(GDof);
F_equiv_local=zeros(GDof,1);
M_Assembly_local=zeros(GDof);
for iElement=1:N_elements
    i_nodes=elementNodes(iElement,:);
    elementDof=[2*(i_nodes(1)-1)+1 2*(i_nodes(2)-1) 2*(i_nodes(2)-1)+1 2*(i_nodes(2)-1)+2];

    L=x_col(i_nodes(2))-x_col(i_nodes(1));

    % Assembly stiffness matrix
    k_local=E*I/L^3*[   12   6*L -12 6*L;
                        6*L 4*L^2 -6*L 2*L^2;
                        -12 -6*L 12 -6*L ;
                        6*L 2*L^2 -6*L 4*L^2];

    K_Assembly_local(elementDof,elementDof)=K_Assembly_local(elementDof,elementDof)+k_local;

    f_eq_1ocal=[ P*L/2
            P*L*L/12
            P*L/2
            -P*L*L/12];

    % Assembly mass matrix
    m_local=rho*A*L/420*[   156 22*L 54 -13*L;
                            22*L 4*L^2 13*L -3*L^2;
                            54 13*L 156 -22*L ;
                            -13*L -3*L^2 -22*L 4*L^2];
    M_Assembly_local(elementDof,elementDof)=M_Assembly_local(elementDof,elementDof)+m_local;

    % equivalent force vector
    F_equiv_local(elementDof)=F_equiv_local(elementDof)+f_eq_1ocal;  
end
