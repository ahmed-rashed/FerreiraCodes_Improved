function K_G_Assembly=formGeometricStiffnessMindlinQ4(GDof,elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness)

% computation of geometric stiffness
% for Mindlin plate element

% K_G_Assembly : geometric matrix
K_G_Assembly=zeros(GDof); 

% Gauss quadrature for bending part
[gaussWeights_reduced,gaussLocations_reduced_cols]=gaussQuadrature(1);
 
N_elements=size(elementNodes,1);
N_nodesPerElement=size(elementNodes,2);
for iElement=1:N_elements       
    i_nodes=elementNodes(iElement,:);           
    elementDof=[i_nodes i_nodes+numberNodes i_nodes+2*numberNodes];    

    % cycle for Gauss point
    for n=1:size(gaussWeights_reduced,1)                      
        xi_reduced_Gauss=gaussLocations_reduced_cols(n,1);
        eta_reduced_Gauss=gaussLocations_reduced_cols(n,2);

        % shape functions and derivatives
        [~,N_diff_xi_eta_cols]=shapeFunctionQ4(xi_reduced_Gauss,eta_reduced_Gauss);
        [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(i_nodes,:),N_diff_xi_eta_cols);

        % geometric matrix
        G_b=zeros(2,3*N_nodesPerElement);
        G_b(1,1:N_nodesPerElement)=N_diff_x_y_cols(:,1)';  
        G_b(2,1:N_nodesPerElement)=N_diff_x_y_cols(:,2)';  
        K_G_Assembly(elementDof,elementDof)=K_G_Assembly(elementDof,elementDof)+G_b'*sigmaMatrix*thickness*G_b*gaussWeights_reduced(n)*det(J_mat);
    end
end

% shear stiffness matrix
for n=1:size(gaussWeights_reduced,1)                      
    xi_reduced_Gauss=gaussLocations_reduced_cols(n,1);
    eta_reduced_Gauss=gaussLocations_reduced_cols(n,2);

    % shape functions and derivatives
    [~,N_diff_xi_eta_cols]=shapeFunctionQ4(xi_reduced_Gauss,eta_reduced_Gauss);
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(i_nodes,:),N_diff_xi_eta_cols);

    % Geometric matrix
    G_s1=zeros(2,3*N_nodesPerElement);
    G_s1(1,N_nodesPerElement+1:2*N_nodesPerElement)=N_diff_x_y_cols(:,1)';  
    G_s1(2,N_nodesPerElement+1:2*N_nodesPerElement)=N_diff_x_y_cols(:,2)';  
    K_G_Assembly(elementDof,elementDof)=K_G_Assembly(elementDof,elementDof)+G_s1'*sigmaMatrix*thickness^3/12*G_s1*gaussWeights_reduced(n)*det(J_mat);

    G_s2=zeros(2,3*N_nodesPerElement);
    G_s2(1,2*N_nodesPerElement+1:3*N_nodesPerElement)=N_diff_x_y_cols(:,1)';  
    G_s2(2,2*N_nodesPerElement+1:3*N_nodesPerElement)=N_diff_x_y_cols(:,2)';  
    K_G_Assembly(elementDof,elementDof)=K_G_Assembly(elementDof,elementDof)+G_s2'*sigmaMatrix*thickness^3/12*G_s2*gaussWeights_reduced(n)*det(J_mat);
end