% Daryl Logan Example 3.5
% antonio ferreira 2008
% Modified by Ahmed Rashed
% Units used are ft lb

clearvars

E_vec=[30e6,30e6,30e6];

A_vec=[2,2,2];

% generation of coordinates and connectivities
nodesCoords=[   0 0;
                0 120;
                120 120;
                120 0];
numberNodes=size(nodesCoords,1);

elementNodes=[1 2;
              1 3;
              1 4];
numberElements=size(elementNodes,1);

% GDof: total number of degrees of freedom
GDof=2*numberNodes;

% Assembly stiffness matrix
K_assembly=formStiffness2Dtruss(GDof,numberElements,elementNodes,nodesCoords,E_vec,A_vec);

% boundary conditions
prescribedDof=3:8;

% force : force vector
F_col=nan(GDof,1);
F_col(1)=0;
F_col(2)=-1e4;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);

% % drawing D_col
% us=1:2:2*numberNodes-1;
% vs=2:2:2*numberNodes;
% figure
% L=nodeCoordinates(2,1)-nodeCoordinates(1,1);
% %L=node(2,1)-node(1,1);
% XX=D_col(us);YY=D_col(vs);
% dispNorm=max(sqrt(XX.^2+YY.^2));
% scaleFact=15000*dispNorm;

% hold on
% drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes,'L2','k.-');
% drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');

% stresses at elements
stress=stresses2Dtruss(numberElements,elementNodes,nodesCoords,D_col,E_vec);
