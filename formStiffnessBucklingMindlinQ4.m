function [K,KG]=...
    formStiffnessBucklingMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,sigmaMatrix,thickness)

% computation of stiffness, mass matrices and force vector 
% for Mindlin plate element
% initial matrices
% K : stiffness matrix
% KG : geometric matrix
K=zeros(GDof);
KG=zeros(GDof); 

% Gauss quadrature for bending part
[W,Q]=gaussQuadrature('reduced');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes];    
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(W,1)                      
    GaussPoint=gaussLocations(q,:);   
    GaussWeight=gaussWeights(q);
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [N,dNdxi]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J0,invJ0,dNdx]=Jacobian(nodeCoordinates(indice,:),dNdxi);
    
% [B] matrix bending
    B_b=zeros(3,3*ndof);
    B_b(1,ndof+1:2*ndof)  = dNdx(:,1)';  
    B_b(2,2*ndof+1:3*ndof)= dNdx(:,2)';
    B_b(3,ndof+1:2*ndof)  = dNdx(:,2)';  
    B_b(3,2*ndof+1:3*ndof)= dNdx(:,1)';
    
% stiffness matrix bending
    K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
        B_b'*C_bending*B_b*W(q)*det(J0);
     
% geometric matrix
    G_b=zeros(2,3*ndof);
    G_b(1,1:ndof)  = dNdx(:,1)';  
    G_b(2,1:ndof)  = dNdx(:,2)';  
    KG(elementDof,elementDof)=KG(elementDof,elementDof)+ ...
        G_b'*sigmaMatrix*thickness*G_b*W(q)*det(J0);
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
  elementDof=[ indice indice+numberNodes indice+2*numberNodes];
  ndof=length(indice);

  % cycle for Gauss point
  for q=1:size(W,1)                        
    pt=Q(q,:);                             
    wt=W(q);         
    xi=pt(1);
    eta=pt(2);
    
% shape functions and derivatives
    %[N,dNdxi]=shapeFunction(xi,eta);
    [N,dNdxi]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [J0,invJ0,dNdx]=Jacobian(nodeCoordinates(indice,:),dNdxi);
    
% [B] matrix shear
    B_s=zeros(2,3*ndof);
    B_s(1,1:ndof)       = dNdx(:,1)';  
    B_s(2,1:ndof)       = dNdx(:,2)';
    B_s(1,ndof+1:2*ndof)  = N;           
    B_s(2,2*ndof+1:3*ndof)= N;

% stiffness matrix shear
    K(elementDof,elementDof)=K(elementDof,elementDof)+B_s'*C_shear  *B_s*W(q)*det(J0);  
    
% Geometric matrix
    G_s1=zeros(2,3*ndof);
    G_s1(1,ndof+1:2*ndof)    = dNdx(:,1)';  
    G_s1(2,ndof+1:2*ndof)    = dNdx(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s1'*sigmaMatrix*thickness^3/12*G_s1*W(q)*det(J0);
    
    G_s2=zeros(2,3*ndof);
    G_s2(1,2*ndof+1:3*ndof)    = dNdx(:,1)';  
    G_s2(2,2*ndof+1:3*ndof)    = dNdx(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s2'*sigmaMatrix*thickness^3/12*G_s2*W(q)*det(J0);

  end  % gauss point
end    % element