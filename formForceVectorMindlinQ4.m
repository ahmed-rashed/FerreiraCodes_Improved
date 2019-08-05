function [force]=formForceVectorMindlinQ4(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,P)

% computation of force vector 
% for Mindlin plate element

% force : force vector
force=zeros(GDof,1);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('reduced');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  indice=elementNodes(e,:);           
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                            
    GaussWeight=gaussWeights(q);                            
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);

% force vector
force(indice)=force(indice)+N_col*P*det(J_mat)*GaussWeight;    
  end  % Gauss point
  
end    % element

