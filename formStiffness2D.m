function [stiffness,mass]=formStiffness2D(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,C,rho,thickness)

% compute stiffness matrix (and mass matrix)
% for plane stress Q4 elements

stiffness=zeros(GDof);
mass=zeros(GDof);

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

for e=1:numberElements
  indice=elementNodes(e,:); 
  elementDof=[ indice indice+numberNodes ];   
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)
    GaussPoint=gaussLocations(q,:);
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta)

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

%  B matrix
    B=zeros(3,2*ndof);
    B(1,1:ndof)       = XYderivatives(:,1)';
    B(2,ndof+1:2*ndof)  = XYderivatives(:,2)';
    B(3,1:ndof)       = XYderivatives(:,2)';
    B(3,ndof+1:2*ndof)  = XYderivatives(:,1)';
    
% stiffness matrix
    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+...
        B'*C*thickness*B*gaussWeights(q)*det(Jacob);
% mass matrix
    mass(indice,indice)=mass(indice,indice)+...
        shapeFunction*shapeFunction'*...
        rho*thickness*gaussWeights(q)*det(Jacob);
    mass(indice+numberNodes,indice+numberNodes)=...
        mass(indice+numberNodes,indice+numberNodes)+...
        shapeFunction*shapeFunction'*...
        rho*thickness*gaussWeights(q)*det(Jacob);

  end
end
