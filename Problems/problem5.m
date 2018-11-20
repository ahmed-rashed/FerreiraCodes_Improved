% Daryl Logan problem 3.32
% antonio ferreira 2008
% Modified by Ahmed Rashed

close all
clearvars

%Units used are N mm

E_vec=70000*ones(1,11);
A_vec=300*ones(1,11);

% generation of coordinates and connectivities
nodesCoords=[   0 0
                0 3000
                3000 0
                3000 3000
                6000 0
                6000 3000];
numberNodes=size(nodesCoords,1);

elementNodes=[  1 2
                1 3
                2 3
                2 4
                1 4
                3 4
                3 6
                4 5
                4 6
                3 5
                5 6];
numberElements=size(elementNodes,1);

% GDof: total number of degrees of freedom
GDof=2*numberNodes;

% Assembly stiffness matrix
K_assembly=formStiffness2Dtruss(GDof,numberElements,elementNodes,nodesCoords,E_vec,A_vec);

% boundary conditions
prescribedDof=[1 2 10];

% force : force vector
F_col=nan(GDof,1);
F_col(4)=-50000;
F_col(8)=-100000;
F_col(12)=-50000;
F_col([3,5,6,7,9,10,11])=0;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);


% % drawing D_col
% us=1:2:2*numberNodes-1;
% vs=2:2:2*numberNodes;
% 
% figure
% L=xx(2)-xx(1);
% %L=node(2,1)-node(1,1);
% XX=D_col(us);YY=D_col(vs);
% dispNorm=max(sqrt(XX.^2+YY.^2));
% scaleFact=2*dispNorm;
% clf
% hold on
% drawingMesh(nodesCoords+scaleFact*[XX YY],...
%     elementNodes,'L2','k.-');
% drawingMesh(nodesCoords,elementNodes,'L2','k.--');
% 
% % output D_col/reactions
% outputDisplacementsReactions(D_col,K_assembly,...
%     GDof,prescribedDof);

% stresses at elements
sigma=stresses2Dtruss(numberElements,elementNodes,nodesCoords,D_col,E_vec);
