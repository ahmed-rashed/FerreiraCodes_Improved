function SrinivasStress(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,qbarra,U,h)

% computes normal and shear stresses forSrinivas case
% note that transverse shear stresses are not corrected

% normal stresses in each layer
  stress_layer1=zeros(numberElements,4,3);
  stress_layer2=zeros(numberElements,4,3);
  stress_layer3=zeros(numberElements,4,3);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature(2);
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % indiceB: element degrees of freedom
  indice=elementNodes(e,:);           
  indiceB=[ indice indice+numberNodes indice+2*numberNodes ...
            indice+3*numberNodes indice+4*numberNodes]; 
  nn=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    pt=gaussLocations(q,:);                            
    wt=gaussWeights(q);                            
    xi=pt(1);
    eta=pt(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);
    
% [B] matrix bending
    B_b=zeros(3,5*nn);
    B_b(1,nn+1:2*nn)        = N_diff_x_y_cols(:,1)';  
    B_b(2,2*nn+1:3*nn)      = N_diff_x_y_cols(:,2)';
    B_b(3,nn+1:2*nn)        = N_diff_x_y_cols(:,2)';  
    B_b(3,2*nn+1:3*nn)      = N_diff_x_y_cols(:,1)';
% [B] matrix membrane
    B_m=zeros(3,5*nn);
    B_m(1,3*nn+1:4*nn)      = N_diff_x_y_cols(:,1)';  
    B_m(2,4*nn+1:5*nn)      = N_diff_x_y_cols(:,2)';
    B_m(3,3*nn+1:4*nn)      = N_diff_x_y_cols(:,2)';  
    B_m(3,4*nn+1:5*nn)      = N_diff_x_y_cols(:,1)';
    
% stresses
    stress_layer1(e,q,:)=...
        2*h/5*qbarra(1:3,1:3,2)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,2)*B_m*U(indiceB);
    stress_layer2(e,q,:)=...
        2*h/5*qbarra(1:3,1:3,3)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,3)*B_m*U(indiceB);
    stress_layer3(e,q,:)=...
        h/2*qbarra(1:3,1:3,3)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,3)*B_m*U(indiceB);

    end  % Gauss point
end    % element

% shear stresses in each layer
% by constitutive equations

    shear_layer1=zeros(numberElements,1,2);
    shear_layer2=zeros(numberElements,1,2);
    shear_layer3=zeros(numberElements,1,2);

% Gauss quadrature for shear part
[gaussWeights,gaussLocations]=gaussQuadrature(1);
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % indiceB: element degrees of freedom
  indice=elementNodes(e,:);           
  indiceB=[ indice indice+numberNodes indice+2*numberNodes ...
            indice+3*numberNodes indice+4*numberNodes]; 
  nn=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    pt=gaussLocations(q,:);                            
    wt=gaussWeights(q);                            
    xi=pt(1);
    eta=pt(2);

% shape functions and derivatives
    [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates(indice,:),N_diff_xi_eta_cols);  
    
% [B] matrix shear
    B_s=zeros(2,5*nn);
    B_s(1,1:nn)       = N_diff_x_y_cols(:,1)';  
    B_s(2,1:nn)       = N_diff_x_y_cols(:,2)';
    B_s(1,nn+1:2*nn)  = N_col;           
    B_s(2,2*nn+1:3*nn)= N_col;

    shear_layer1(e,q,:)=qbarra(4:5,4:5,1)*B_s*U(indiceB);
    shear_layer2(e,q,:)=qbarra(4:5,4:5,2)*B_s*U(indiceB);
  
  end  % gauss point
end    % element

% normalized stresses, look for table in the book
format
[ abs(min(stress_layer3(:,3,1))),...
  abs(min(stress_layer2(:,3,1))), ...
  abs(min(stress_layer1(:,3,1))),...
  max(shear_layer2(:,:,1)),...
  max(shear_layer1(:,:,1))]

