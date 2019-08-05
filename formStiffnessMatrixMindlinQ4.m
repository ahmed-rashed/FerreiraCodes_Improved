%................................................................

function [K]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I)

% computation of stiffness matrix
% for Mindlin plate element

% K : stiffness matrix

K=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes];    
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
    
% stiffness matrix bending
    K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
        B_b'*C_bending*B_b*gaussWeights(q)*det(J_mat);
    end  % Gauss point
end    % element

% shear stiffness matrix

% Gauss quadrature for shear part
[gaussWeights,gaussLocations]=gaussQuadrature('reduced');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes];    
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
% [B] matrix shear
    B_s=zeros(2,3*ndof);
    B_s(1,1:ndof)       = N_diff_x_y_cols(:,1)';  
    B_s(2,1:ndof)       = N_diff_x_y_cols(:,2)';
    B_s(1,ndof+1:2*ndof)  = N_col;           
    B_s(2,2*ndof+1:3*ndof)= N_col;

% stiffness matrix shear
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
        B_s'*C_shear  *B_s*gaussWeights(q)*det(J_mat);  
  end  % gauss point
end    % element
