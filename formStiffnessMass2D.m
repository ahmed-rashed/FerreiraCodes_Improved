% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira
% plane stress Q4 elements

function [K_Assembly,M_Assembly]=formStiffnessMass2D(GDof,elementNodes,nodeCoordinates,D,rho,h)
K_Assembly=zeros(GDof,GDof);
M_Assembly=zeros(GDof,GDof);

% 2 by 2 quadrature
[gaussWeights,gaussLocations_cols]=gaussQuadrature('complete');

N_elements=size(elementNodes,1);
N_nodesPerElement=size(elementNodes,2);
elementDof=nan(1,2*N_nodesPerElement);
B=zeros(3,2*N_nodesPerElement);
for iElement=1:N_elements
    i_nodes=elementNodes(iElement,:); 
    elementDof(1:2:end)=2*i_nodes-1;
    elementDof(2:2:end)=2*i_nodes;

    N_mat=zeros(2,8);
    k_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    m_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    % cycle for Gauss point
    for n=1:size(gaussWeights,1)
        xi_Gauss=gaussLocations_cols(n,1);
        eta_Gauss=gaussLocations_cols(n,2);

        [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi_Gauss,eta_Gauss);
        N_mat(1,1:2:7)=N_col.';
        N_mat(1,2:2:8)=N_col.';
        [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(i_nodes,:),N_diff_xi_eta_cols);

        B(1,1:2:end)=N_diff_x_y_cols(:,1).';
        B(2,2:2:end)=N_diff_x_y_cols(:,2).';
        B(3,1:2:end)=N_diff_x_y_cols(:,2).';
        B(3,2:2:end)=N_diff_x_y_cols(:,1).';

        k_element=k_element+B.'*D*h*B*gaussWeights(n)*det(J_mat);
        m_element=m_element+h*rho*(N_mat.'*N_mat)*gaussWeights(n)*det(J_mat);
    end
    K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+k_element;
    M_Assembly(elementDof,elementDof)=M_Assembly(elementDof,elementDof)+m_element;
end