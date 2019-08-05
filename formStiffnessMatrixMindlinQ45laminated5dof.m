%................................................................

function [K]=...
    formStiffnessMatrixMindlinQ45laminated5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,AMatrix,...
    BMatrix,DMatrix,SMatrix)

% computation of stiffness matrix
% for Mindlin plate element

% K : stiffness matrix

K=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[ indice indice+numberNodes indice+2*numberNodes ...
            indice+3*numberNodes indice+4*numberNodes]; 
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);
    
% [B] matrix bending
    B_b=zeros(3,5*ndof);
    B_b(1,ndof+1:2*ndof)        = N_diff_x_y_cols(:,1)';  
    B_b(2,2*ndof+1:3*ndof)      = N_diff_x_y_cols(:,2)';
    B_b(3,ndof+1:2*ndof)        = N_diff_x_y_cols(:,2)';  
    B_b(3,2*ndof+1:3*ndof)      = N_diff_x_y_cols(:,1)';
% [B] matrix membrane
    B_m=zeros(3,5*ndof);
    B_m(1,3*ndof+1:4*ndof)      = N_diff_x_y_cols(:,1)';  
    B_m(2,4*ndof+1:5*ndof)      = N_diff_x_y_cols(:,2)';
    B_m(3,3*ndof+1:4*ndof)      = N_diff_x_y_cols(:,2)';  
    B_m(3,4*ndof+1:5*ndof)      = N_diff_x_y_cols(:,1)';
    
% stiffness matrix 

% ... bending-bending
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
                       B_b'*DMatrix*B_b*gaussWeights(q)*det(J_mat);
% ... membrane-membrane                  
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
                       B_m'*AMatrix*B_m*gaussWeights(q)*det(J_mat);
% ... membrane-bending                  
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
                       B_m'*BMatrix*B_b*gaussWeights(q)*det(J_mat);
% ... bending-membrane                  
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
                       B_b'*BMatrix*B_m*gaussWeights(q)*det(J_mat);

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
  elementDof=[ indice indice+numberNodes indice+2*numberNodes ...
            indice+3*numberNodes indice+4*numberNodes]; 
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);  
    
% [B] matrix shear
    B_s=zeros(2,5*ndof);
    B_s(1,1:ndof)       = N_diff_x_y_cols(:,1)';  
    B_s(2,1:ndof)       = N_diff_x_y_cols(:,2)';
    B_s(1,ndof+1:2*ndof)  = N_col;           
    B_s(2,2*ndof+1:3*ndof)= N_col;

% stiffness matrix shear
    K(elementDof,elementDof)=K(elementDof,elementDof)+...
        B_s'*SMatrix  *B_s*gaussWeights(q)*det(J_mat);  
  end  % gauss point
end    % element
