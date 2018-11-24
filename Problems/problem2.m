% antonio ferreira 2008
% Modified by Ahmed Rashed

clearvars

E=30e6;
A=1;

% generation of coordinates and connectivities
L=90;
N_nodes=4;
x_col=linspace(0,L,N_nodes).';

N_elements=N_nodes-1;
elementNodes=[(1:N_elements).',(2:N_elements+1).'];

% GDof: global number of degrees of freedom
GDof=2*N_nodes;

% computation of the system stiffness matrix
K_Assembly=zeros(GDof,GDof);
for iElement=1:N_elements
    elementDof=elementNodes(iElement,:);
    N_elementDof=length(elementDof);    
    L_element=x_col(elementDof(2))-x_col(elementDof(1));
    detJacobian=L_element/2;
    invJacobian=1/detJacobian;

    % central Gauss point (xi=0, weight W=2)
    [shape,naturalDerivatives]=shapeFunctionL2(0);
    Xderivatives=naturalDerivatives*invJacobian;

    % B matrix     
    B_row=Xderivatives(:).'; 
    K_Assembly(elementDof,elementDof)=K_Assembly(elementDof,elementDof)+B_row'*B_row*2*detJacobian*E*A;
end 

D_col=zeros(numberNodes,1);
F_col=zeros(numberNodes,1);
K_Assembly=zeros(numberNodes,numberNodes); 
% applied load at node 2
F_col(2)=3000.0;

% boundary conditions and solution
% prescribed dofs
fixedDof=find(x_col==min(x_col(:)) | x_col==max(x_col(:)))'; 
prescribedDof=fixedDof
% free Dof : activeDof
activeDof=setdiff([1:numberNodes]',prescribedDof);

% solution
GDof=numberNodes;
D_col=solution(GDof,prescribedDof,K_Assembly,F_col);

% output D_col/reactions
outputDisplacementsReactions(D_col,K_Assembly,numberNodes,prescribedDof)
