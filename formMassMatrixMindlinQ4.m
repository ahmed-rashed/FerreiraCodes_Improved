function [mass]=formMassMatrixMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,thickness,rho,I)

% computation of  mass matrix 
% for Mindlin plate element

% mass : mass matrix 
mass=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
 
% cycle for element
for e=1:numberElements
  % indice : nodal condofectivities for each element
  indice=elementNodes(e,:);
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)
    GaussPoint=gaussLocations(q,:);
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);
    
% [B] matrix bending
    B_b=zeros(3,3*ndof);
    B_b(1,ndof+1:2*ndof)  = N_diff_x_y_cols(:,1)';  
    B_b(2,2*ndof+1:3*ndof)= N_diff_x_y_cols(:,2)';
    B_b(3,ndof+1:2*ndof)  = N_diff_x_y_cols(:,2)';  
    B_b(3,2*ndof+1:3*ndof)= N_diff_x_y_cols(:,1)';
    
% mass matrix
    
    mass(indice,indice)=mass(indice,indice)+N_col*N_col'*thickness*rho*gaussWeights(q)*det(J_mat);
    mass(indice+numberNodes,indice+numberNodes)=mass(indice+numberNodes,indice+numberNodes)+N_col*N_col'*I*rho*gaussWeights(q)*det(J_mat);
    mass(indice+2*numberNodes,indice+2*numberNodes)=mass(indice+2*numberNodes,indice+2*numberNodes)+N_col*N_col'*I*rho*gaussWeights(q)*det(J_mat);
   
  end  % Gauss point
end    % element

