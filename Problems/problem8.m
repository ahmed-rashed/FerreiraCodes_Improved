% 3D truss
% Daryl Logan
% antonio ferreira 2008
% Modified by Ahmed Rashed

close all
clearvars


E_vec=210000*ones(1,4); 
A_vec=100*ones(1,4);

% generation of coordinates and connectivities
nodesCoords=[   4000 4000 3000;
                0 4000 0;
                0 4000 6000;
                4000 0    3000;
                8000 -1000 1000];
numberNodes=size(nodesCoords,1);

elementNodes=[  1 2
                1 3
                1 4
                1 5];
numberElements=size(elementNodes,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes; 

% Assembly stiffness matrix
K_assembly=formStiffness3Dtruss(GDof,numberElements,elementNodes,nodesCoords,E_vec,A_vec);



% boundary conditions
prescribedDof=[4:15];

%force vector
F_col=nan(GDof,1);
F_col([1,3])=0;
F_col(2)=-10000;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);

% % output D_col/reactions
% outputDisplacementsReactions(D_col,K_assembly,GDof,prescribedDof);

% stresses at elements
stresses3Dtruss(numberElements,elementNodes,nodesCoords,D_col,E_vec);
