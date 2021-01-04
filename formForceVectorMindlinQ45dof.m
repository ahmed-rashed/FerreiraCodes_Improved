function force=formForceVectorMindlinQ45dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,P)

% computation of force vector 
% for Mindlin plate element

% force : force vector
force=zeros(GDof,1);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature(1);
 
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
    [~,N_diff_xi_eta_rows]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    J_mat=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_rows);
    
% force vector
    force(indice)=force(indice)+shapeFunction*P*det(J_mat)*GaussWeight;
  end  % Gauss point
  
end    % element

