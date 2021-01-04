% Modified by Ahmed Rashed

function [K_Assembly,M_Assembly,F_eq,K_G_Assembly]=formMatricesMindlinQ4(GDof,elementNodes,nodeCoordinates,h,D_s,D_b,p_zeta)

% 2x2 Gauss quadrature for bending part
[gaussWeights,gaussLocations_cols]=gaussQuadrature(2);

% 1x1 Gauss quadrature for shear part
[gaussWeights_reduced,gaussLocations_reduced_cols]=gaussQuadrature(1);

N_elements=size(elementNodes,1);
N_nodesPerElement=size(elementNodes,2);
elementDof=nan(1,3*N_nodesPerElement);
B_b=zeros(3,3*N_nodesPerElement);
K_Assembly=zeros(GDof,GDof);
M_Assembly=zeros(GDof,GDof);
F_eq=zeros(GDof,1);
r_diag=diag([h,h^3/12,h^3/12]);
for iElement=1:N_elements
    i_nodes=elementNodes(iElement,:);           
    elementDof(1:3:end)=3*i_nodes-2;
    elementDof(2:3:end)=3*i_nodes-1;
    elementDof(3:3:end)=3*i_nodes;

    N_mat=zeros(3,3*N_nodesPerElement);
    k_b_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    m_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    % bending stiffness matrix
    for n=1:size(gaussWeights,1)    % cycle for Gauss point
        xi_Gauss=gaussLocations_cols(n,1);
        eta_Gauss=gaussLocations_cols(n,2);

        [N_col,N_diff_xi_eta_rows]=shapeFunctionQ4(xi_Gauss,eta_Gauss);
        N_mat(1,1:3:end)=N_col.';
        N_mat(2,2:3:end)=N_col.';
        N_mat(3,3:3:end)=N_col.';
        [J_mat,N_diff_x_y_rows]=Jacobian(nodeCoordinates(i_nodes,:),N_diff_xi_eta_rows);

        % [B_b] matrix bending
        B_b(1,3:3:end)=N_diff_x_y_rows(1,:);
        B_b(2,2:3:end)=-N_diff_x_y_rows(2,:);
        B_b(3,2:3:end)=-N_diff_x_y_rows(1,:);
        B_b(3,3:3:end)=N_diff_x_y_rows(2,:);

        k_b_element=k_b_element+h^3/12*B_b.'*D_b*B_b*gaussWeights(n)*det(J_mat);
        m_element=m_element+rho*(N_mat.'*r_diag*N_mat)*gaussWeights(n)*det(J_mat);
    end
    M_Assembly(elementDof,elementDof)=M_Assembly(elementDof,elementDof)+m_element;
    
    % shear stiffness matrix
    k_s_element=zeros(2*N_nodesPerElement,2*N_nodesPerElement);
    f_eq_element=zeros(N_nodesPerElement,1);
    for n=1:size(gaussWeights_reduced,1)    % cycle for Gauss point
        xi_reduced_Gauss=gaussLocations_reduced_cols(n,1);
        eta_reduced_Gauss=gaussLocations_reduced_cols(n,2);

        % shape functions and derivatives
        [N_col,N_diff_xi_eta_rows]=shapeFunctionQ4(xi_reduced_Gauss,eta_reduced_Gauss);

        % Jacobian matrix, inverse of Jacobian, 
        % derivatives w.r.t. x,y    
        [J_mat,N_diff_x_y_rows]=Jacobian(nodeCoordinates(i_nodes,:),N_diff_xi_eta_rows);    
        % [B] matrix shear
        B_s=zeros(2,3*N_nodesPerElement);
        B_s(1,1:N_nodesPerElement)       = N_diff_x_y_rows(1,:);
        B_s(2,1:N_nodesPerElement)       = N_diff_x_y_rows(2,:);
        B_s(1,N_nodesPerElement+1:2*N_nodesPerElement)  = N_col;
        B_s(2,2*N_nodesPerElement+1:3*N_nodesPerElement)= N_col;

        k_s_element=k_s_element+B_s'*D_s*B_s*det(J_mat)*gaussWeights_reduced(n);
        f_eq_element=f_eq_element+N_col*p_zeta*det(J_mat)*gaussWeights_reduced(n);
    end
    K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+k_b_element+k_s_element;
    F_eq(elementDof(1:3:end))=F_eq(elementDof(1:3:end))+f_eq_element;
end