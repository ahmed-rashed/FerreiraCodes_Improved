% Daryl Logan Example 3.5
% Units used are lb in
% antonio ferreira 2008
% Modified by Ahmed Rashed

close all
clearvars

E_vec=1.2e6*ones(1,3); 
A_vec=[0.302;0.729;0.187]; % area for various sections

% generation of coordinates and connectivities
nodesCoords=[   72 0 0
                0 36 0
                0 36 72
                0 0 -48];
numberNodes=size(nodesCoords,1);

elementNodes=[  1 2
                1 3
                1 4]; 
numberElements=size(elementNodes,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% Assembly stiffness matrix
K_assembly=formStiffness3Dtruss(GDof,numberElements,elementNodes,nodesCoords,E_vec,A_vec);

% boundary conditions
prescribedDof=[2 4:12];

% force : force vector
F_col=nan(GDof,1);
F_col([1,2])=0;
F_col(3)=-1e3;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);

% % output D_col/reactions
% outputDisplacementsReactions(D_col,K_assembly,...
%     GDof,prescribedDof)

% stresses at elements
stresses3Dtruss(numberElements,elementNodes,nodesCoords,D_col,E_vec);

