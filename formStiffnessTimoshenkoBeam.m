function [stiffness,force]=...
    formStiffnessTimoshenkoBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,C,P)

% computation of stiffness matrix and force vector 
% for Timoshenko beam element
stiffness=zeros(GDof);
force=zeros(GDof,1);

    gaussLocations=[0.577350269189626;-0.577350269189626];                            
    gaussWeights=ones(2,1);
    
% bending contribution for stiffness matrix
for e=1:numberElements
    indice=elementNodes(e,:);
    elementDof=[ indice indice+numberNodes];
    ndof=length(indice);    
    length_element=xx(indice(2))-xx(indice(1));
    detJacobian=length_element/2;invJacobian=1/detJacobian;
  for q=1:size(gaussWeights,1) ; 
     pt=gaussLocations(q,:); 
     [shape,naturalDerivatives]=shapeFunctionL2(pt(1));
     Xderivatives=naturalDerivatives*invJacobian;
% B matrix
B=zeros(2,2*ndof);  
B(1,ndof+1:2*ndof)  = Xderivatives(:)'; 
% K
stiffness(elementDof,elementDof)=...
    stiffness(elementDof,elementDof)+...
    B'*B*gaussWeights(q)*detJacobian*C(1,1);
force(indice)=force(indice)+shape*P*detJacobian*gaussWeights(q);    
  end  
end   
 
% shear component
    gaussLocations=[0.];                            
    gaussWeights=[1.];

for e=1:numberElements  
    indice=elementNodes(e,:); 
    elementDof=[ indice indice+numberNodes];
ndof=length(indice);  
length_element=xx(indice(2))-xx(indice(1));
detJ0=length_element/2;invJ0=1/detJ0;
  for q=1:size(gaussWeights,1) ; 
     pt=gaussLocations(q,:);
     [shape,naturalDerivatives]=shapeFunctionL2(pt(1));
     Xderivatives=naturalDerivatives*invJacobian;
% B     
    B=zeros(2,2*ndof);    
    B(2,1:ndof)       = Xderivatives(:)';  
    B(2,ndof+1:2*ndof)  = shape;
% K
    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+...
        B'*B*gaussWeights(q)*detJacobian*C(2,2); 
  end  
end 