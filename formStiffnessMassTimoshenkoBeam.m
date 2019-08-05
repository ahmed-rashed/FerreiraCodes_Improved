function [K_Assembly,F_equiv_local,M_Assembly]=formStiffnessMassTimoshenkoBeam(GDof,elementNodes,x_col,C,P,rho,I,h)

N_elements=size(elementNodes,1);
N_Nodes=length(x_col);

K_Assembly=zeros(GDof,GDof);
M_Assembly=zeros(GDof,GDof);
F_equiv_local=zeros(GDof,1);

% stiffness matrix
gaussLocations=[0.577350269189626
                -0.577350269189626];                            
gaussWeights=ones(2,1);
    
% bending contribution for stiffness matrix
for iElement=1:N_elements
    i_nodes=elementNodes(iElement,:);
    elementDof=[i_nodes i_nodes+N_Nodes];
    indiceMass=i_nodes+N_Nodes;
    ndof=length(i_nodes);    
    L=x_col(i_nodes(2))-x_col(i_nodes(1));
    detJacobian=L/2;
    invJacobian=1/detJacobian;
    for q=1:size(gaussWeights,1)
        pt=gaussLocations(q,:); 
        [shape,naturalDerivatives]=shapeFunctionL2(pt(1));
        Xderivatives=naturalDerivatives*invJacobian;
        
        % B matrix
        B=zeros(2,2*ndof);  
        B(1,ndof+1:2*ndof)=Xderivatives(:)'; 
        
        % K_Assembly
        K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+B'*B*gaussWeights(q)*detJacobian*C(1,1);
        F_equiv_local(i_nodes)=F_equiv_local(i_nodes)+shape*P*detJacobian*gaussWeights(q); 

        M_Assembly(indiceMass,indiceMass)=M_Assembly(indiceMass,indiceMass) +shape*shape'*gaussWeights(q)*I*rho*detJacobian;
        M_Assembly(i_nodes,i_nodes)=M_Assembly(i_nodes,i_nodes)             +shape*shape'*gaussWeights(q)*h*rho*detJacobian;
    end  
end

% shear contribution for stiffness matrix
gaussLocations=[0.];                            
gaussWeights=[2.];
for iElement=1:N_elements  
    i_nodes=elementNodes(iElement,:); 
    elementDof=[ i_nodes i_nodes+N_Nodes];
    ndof=length(i_nodes);  
    L=x_col(i_nodes(2))-x_col(i_nodes(1));
    detJacobian=L/2;
    invJacobian=1/detJacobian;
    for q=1:size(gaussWeights,1)
        pt=gaussLocations(q,:);
        [shape,naturalDerivatives]=shapeFunctionL2(pt(1));
        Xderivatives=naturalDerivatives*invJacobian;
        
        % B     
        B=zeros(2,2*ndof);    
        B(2,1:ndof)       = Xderivatives(:)';  
        B(2,ndof+1:2*ndof)  = shape;
        
        % K_Assembly
        K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+B'*B*gaussWeights(q)*detJacobian*C(2,2); 
    end  
end 