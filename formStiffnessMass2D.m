% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

function [K_Assembly,M_Assembly]=formStiffnessMass2D(GDof,elementNodes,nodeCoordinates,C,rho,h)

% compute stiffness matrix (and mass matrix) for plane stress Q4 elements

K_Assembly=zeros(GDof,GDof);
M_Assembly=zeros(GDof,GDof);

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

N_elements=size(elementNodes,1);
N_nodesPerElement=size(elementNodes,2);
elementDof=nan(1,2*N_nodesPerElement);
B=zeros(3,2*N_nodesPerElement);
for iElement=1:N_elements
    i_nodes=elementNodes(iElement,:); 
    elementDof(1:2:end)=2*i_nodes-1;
    elementDof(2:2:end)=2*i_nodes;

    k_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    m_element_reduced=zeros(N_nodesPerElement,N_nodesPerElement);
    % cycle for Gauss point
    for q=1:size(gaussWeights,1)
        xi=gaussLocations(q,1);
        eta=gaussLocations(q,2);

        [N_reduced_col,N_reduced_xi_eta_cols]=shapeFunctionQ4(xi,eta);
        [J_mat,N_reduced_x_y_cols]=Jacobian(nodeCoordinates(i_nodes,:),N_reduced_xi_eta_cols);

        B(1,1:2:end)=N_reduced_x_y_cols(:,1).';
        B(2,2:2:end)=N_reduced_x_y_cols(:,2).';
        B(3,1:2:end)=N_reduced_x_y_cols(:,2).';
        B(3,2:2:end)=N_reduced_x_y_cols(:,1).';

        k_element=k_element+B.'*C*h*B*gaussWeights(q)*det(J_mat);
        m_element_reduced=m_element_reduced+N_reduced_col*N_reduced_col.'*rho*h*gaussWeights(q)*det(J_mat);
    end
    K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+k_element;
    M_Assembly(elementDof(1:2:end),elementDof(1:2:end))=M_Assembly(elementDof(1:2:end),elementDof(1:2:end))+m_element_reduced;
    M_Assembly(elementDof(2:2:end),elementDof(2:2:end))=M_Assembly(elementDof(2:2:end),elementDof(2:2:end))+m_element_reduced;
end