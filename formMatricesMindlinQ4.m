function [K,KG,mass,force]=...
    formMatricesMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,sigmaMatrix,thickness,rho,I,P)

% computation of stiffness, mass matrices and force vector 
% for Mindlin plate element
% initial matrices
% K : stiffness matrix
% KG : geometric matrix
% mass : mass matrix
% force : force vector
K=zeros(GDof);
KG=zeros(GDof); 
mass=zeros(GDof);
force=zeros(GDof,1);

% Gauss quadrature for bending part
[W,Q]=gaussQuadrature('reduced');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[indice indice+numberNodes indice+2*numberNodes];    
  nn=length(indice);
  
  % cycle for Gauss point
  for q=1:size(W,1)                      
    pt=Q(q,:);                            
    wt=W(q);                            
    xi=pt(1);
    eta=pt(2);

% shape functions and derivatives
    [N,dNdxi]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,dNdx]=Jacobian(nodeCoordinates(indice,:),dNdxi);
    
% [B] matrix bending
    B_b=zeros(3,3*nn);
    B_b(1,nn+1:2*nn)  = dNdx(:,1)';  
    B_b(2,2*nn+1:3*nn)= dNdx(:,2)';
    B_b(3,nn+1:2*nn)  = dNdx(:,2)';  
    B_b(3,2*nn+1:3*nn)= dNdx(:,1)';
    
% stiffness matrix bending
    K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
        B_b'*C_bending*B_b*W(q)*det(Jacob);
     
% geometric matrix
    G_b=zeros(2,3*nn);
    G_b(1,1:nn)  = dNdx(:,1)';  
    G_b(2,1:nn)  = dNdx(:,2)';  
    KG(elementDof,elementDof)=KG(elementDof,elementDof)+ ...
        G_b'*sigmaMatrix*thickness*G_b*W(q)*det(Jacob);
    
        mass(indice,indice)=mass(indice,indice)+...
        N*N'*thickness*rho*W(q)*det(Jacob);
    mass(indice+numberNodes,indice+numberNodes)=...
        mass(indice+numberNodes,indice+numberNodes)+...
        N*N'*I*rho*W(q)*det(Jacob);
    mass(indice+2*numberNodes,indice+2*numberNodes)=...
        mass(indice+2*numberNodes,indice+2*numberNodes)+...
        N*N'*I*rho*W(q)*det(Jacob);
% force vector
    force(indice)=force(indice)+N*P*det(Jacob)*wt;    
  end  % Gauss point
  
end    % element

% shear stiffness matrix

% Gauss quadrature for shear part
[W,Q]=gaussQuadrature('reduced');

% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % elementDof: element degrees of freedom
  indice=elementNodes(e,:);           
  elementDof=[ indice indice+numberNodes indice+2*numberNodes];
  nn=length(indice);

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
    [Jacob,invJacobian,dNdx]=Jacobian(nodeCoordinates(indice,:),dNdxi);
    
% [B] matrix shear
    B_s=zeros(2,3*nn);
    B_s(1,1:nn)       = dNdx(:,1)';  
    B_s(2,1:nn)       = dNdx(:,2)';
    B_s(1,nn+1:2*nn)  = N;           
    B_s(2,2*nn+1:3*nn)= N;

% stiffness matrix shear
    K(elementDof,elementDof)=K(elementDof,elementDof)+B_s'*C_shear  *B_s*W(q)*det(Jacob);  
    
% Geometric matrix
    G_s1=zeros(2,3*nn);
    G_s1(1,nn+1:2*nn)    = dNdx(:,1)';  
    G_s1(2,nn+1:2*nn)    = dNdx(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s1'*sigmaMatrix*thickness^3/12*G_s1*W(q)*det(Jacob);
    
    G_s2=zeros(2,3*nn);
    G_s2(1,2*nn+1:3*nn)    = dNdx(:,1)';  
    G_s2(2,2*nn+1:3*nn)    = dNdx(:,2)';  
    KG(elementDof,elementDof)  =KG(elementDof,elementDof)+ ...
        G_s2'*sigmaMatrix*thickness^3/12*G_s2*W(q)*det(Jacob);

  end  % gauss point
end    % element