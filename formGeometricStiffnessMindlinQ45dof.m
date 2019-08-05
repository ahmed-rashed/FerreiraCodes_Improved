%................................................................

function [KG]=...
    formGeometricStiffnessMindlinQ45dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness)

% computation of geometric stiffness
% for Mindlin plate element

% KG : geometric matrix
KG=zeros(GDof); 

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('reduced');

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
    [~,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);

% geometric matrix
    G_b=zeros(2,5*ndof);
    G_b(1,1:ndof)  = N_diff_x_y_cols(:,1)';  
    G_b(2,1:ndof)  = N_diff_x_y_cols(:,2)';  
    KG(elementDof,elementDof)=KG(elementDof,elementDof)+ ...
        G_b'*sigmaMatrix*thickness*G_b*gaussWeights(q)*det(J_mat);
    
   end  % Gauss point
  
end    % element

% shear stiffness matrix

% Gauss quadrature for shear part
[W,Q]=gaussQuadrature('reduced');

% cycle for element

for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[ indice indice+numberNodes indice+2*numberNodes ...
            indice+3*numberNodes indice+4*numberNodes]; 
  ndof=length(indice);

  
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [~,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);

        
% Geometric matrix
    G_s1=zeros(2,5*ndof);
    G_s1(1,ndof+1:2*ndof)    = N_diff_x_y_cols(:,1)';  
    G_s1(2,ndof+1:2*ndof)    = N_diff_x_y_cols(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s1'*sigmaMatrix*thickness^3/12*G_s1*gaussWeights(q)*det(J_mat);
    
    G_s2=zeros(2,5*ndof);
    G_s2(1,2*ndof+1:3*ndof)    = N_diff_x_y_cols(:,1)';  
    G_s2(2,2*ndof+1:3*ndof)    = N_diff_x_y_cols(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s2'*sigmaMatrix*thickness^3/12*G_s2*gaussWeights(q)*det(J_mat);

  end  % gauss point
end    % element
