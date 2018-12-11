function [K_Assembly,M_Assembly]=formStiffness2D(GDof,numberElements,elementNodes,nodeCoordinates,C,rho,h)

% compute stiffness matrix (and mass matrix) for plane stress Q4 elements

K_Assembly=zeros(GDof,GDof);
M_Assembly=zeros(GDof,GDof);

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

N_nodesElement=size(elementNodes,2);
for iElement=1:numberElements
    i_nodes=elementNodes(iElement,:); 
    elementDof=[i_nodes i_nodes+N_nodesElement];

    % cycle for Gauss point
    for q=1:size(gaussWeights,1)
        xi=gaussLocations(q,1);
        eta=gaussLocations(q,2);

        % shape functions and derivatives
        [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);

        % Jacobian matrix, inverse of Jacobian, derivatives w.r.t. x,y    
        [Jacob,invJacobian,XYderivatives]=Jacobian(nodeCoordinates(i_nodes,:),naturalDerivatives);

        %  B matrix
        B=zeros(3,2*N_nodesElement);
        B(1,1:N_nodesElement)           = XYderivatives(:,1).';
        B(2,N_nodesElement+1:2*N_nodesElement) = XYderivatives(:,2).';
        B(3,1:N_nodesElement)           = XYderivatives(:,2).';
        B(3,N_nodesElement+1:2*N_nodesElement) = XYderivatives(:,1).';

        K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+B.'*C*h*B*gaussWeights(q)*det(Jacob);

        M_Assembly(i_nodes,i_nodes)=M_Assembly(i_nodes,i_nodes)+shapeFunction*shapeFunction.'*rho*h*gaussWeights(q)*det(Jacob);
        M_Assembly(i_nodes+N_nodesElement,i_nodes+N_nodesElement)=M_Assembly(i_nodes+N_nodesElement,i_nodes+N_nodesElement)+shapeFunction*shapeFunction.'*rho*h*gaussWeights(q)*det(Jacob);
    end
end