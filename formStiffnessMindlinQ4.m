%................................................................

function [stiffness,force]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,C_bending)

% computation of stiffness matrix and force vector 
% for Mindlin plate element
stiffness=zeros(GDof,GDof);

% 2 x 2 Gauss quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % indiceB: element degrees of freedom
  indice=elementNodes(e,:);           
  indiceB=[indice indice+numberNodes indice+2*numberNodes];    
  nn=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    pt=gaussLocations(q,:);                            
    wt=gaussWeights(q);                            
    xi=pt(1);
    eta=pt(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

%  B matrix
    B_b=zeros(3,3*nn);
    B_b(1,nn+1:2*nn)  = XYderivatives(:,1)';  
    B_b(2,2*nn+1:3*nn)= XYderivatives(:,2)';
    B_b(3,nn+1:2*nn)  = XYderivatives(:,2)';  
    B_b(3,2*nn+1:3*nn)= XYderivatives(:,1)';
% stiffness matrix
    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+...
        B_b'*C_bending*B_b*W(q)*det(jacobian);    
  end  
end    
%  stiffness matrix (shear contribution)
%  1 point Gauss quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('reduced');
 
% cycle for element
% cycle for element
for e=1:numberElements       
  % indice : nodal connectivities for each element
  % indiceB: element degrees of freedom
  indice=elementNodes(e,:);           
  indiceB=[indice indice+numberNodes indice+2*numberNodes];    
  nn=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    pt=gaussLocations(q,:);                            
    wt=gaussWeights(q);                            
    xi=pt(1);
    eta=pt(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

%  B matrix
    B_s=zeros(2,3*nn);
    B_s(1,1:nn)       = XYderivatives(:,1)';  
    B_s(2,1:nn)       = XYderivatives(:,2)';
    B_s(1,nn+1:2*nn)  = shapeFunction;
    B_s(2,2*nn+1:3*nn)= shapeFunction;   
% stiffness matrix
    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+...
        B_s'*C_shear  *B_s*W(q)*det(jacobian);    
  end  
end    
