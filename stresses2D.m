function stresses=stresses2D(numberElements,elementNodes,numberNodes,nodeCoordinates,D_col,UX,UY,C,scaleFactor)

% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');

% stresses at nodes
stresses=zeros(numberElements,size(elementNodes,2),3);
  
for iElement=1:numberElements                           
    i_nodes=elementNodes(iElement,:); 
    elementDof=[ i_nodes i_nodes+numberNodes ];   
    nn=length(i_nodes);
    for q=1:size(gaussWeights,1)                        
        xi=gaussLocations(q1);
        eta=gaussLocations(q2);
       
        [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);

        [Jacob,XYderivatives]=Jacobian(nodeCoordinates(i_nodes,:),naturalDerivatives);

        %  B matrix
        B=zeros(3,2*nn);
        B(1,1:nn)       = XYderivatives(:,1)';
        B(2,nn+1:2*nn)  = XYderivatives(:,2)';
        B(3,1:nn)       = XYderivatives(:,2)';
        B(3,nn+1:2*nn)  = XYderivatives(:,1)';

        % element deformation 
        strain=B*D_col(elementDof);
        stresses(iElement,q,:)=C*strain;    
    end                               
end   

% % drawing stress fields
% % on top of the deformed shape
% figure
% drawingField(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4',stress(:,:,1));%sigma XX
% hold on
% drawingMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4','k-');
% drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
% colorbar
% title('Sigma XX stress (on deformed shape)')
% axis off
