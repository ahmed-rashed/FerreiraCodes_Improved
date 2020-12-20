function stresses=stresses2D(elementNodes,N_nodes,nodeCoordinates,D_col,UX,UY,C,scaleFactor)

N_elements=size(elementNodes,1);

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature(2);

% stresses at nodes
N_nodesPerElement=size(elementNodes,2);
stresses=zeros(N_elements,3,N_nodesPerElement);

elementDof=nan(1,2*N_nodesPerElement);
B=zeros(3,2*N_nodesPerElement);
for iElement=1:N_elements                           
    i_nodes=elementNodes(iElement,:); 
    elementDof(1:2:end)=2*i_nodes-1;
    elementDof(2:2:end)=2*i_nodes;
    for q=1:size(gaussWeights,1)
        xi=gaussLocations(q,1);
        eta=gaussLocations(q,2);
       
        [~,N_reduced_xi_eta_cols]=shapeFunctionQ4(xi,eta);
        [~,N_diff_x_y_cols]=Jacobian(nodeCoordinates(i_nodes,:),N_reduced_xi_eta_cols);

        B(1,1:2:end)=N_diff_x_y_cols(:,1).';
        B(2,2:2:end)=N_diff_x_y_cols(:,2).';
        B(3,1:2:end)=N_diff_x_y_cols(:,2).';
        B(3,2:2:end)=N_diff_x_y_cols(:,1).';

        stresses(iElement,:,q)=C*B*D_col(elementDof);
    end                               
end   

% % drawing stress fields on top of the deformed shape
% figure
% drawingField(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4',stresses(:,:,1));%sigma XX
% hold on
% drawingMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4','k-');
% drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
% colorbar
% title('Sigma XX stress (on deformed shape)')
% axis off
