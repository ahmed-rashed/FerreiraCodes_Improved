% antonio ferreira 2008
%Modified by Ahmed Rashed
%This corrects the strange node ordering of Ferreira's book

clearvars

E_vec=210000*[1,1,1];
A_vec=100*[1,1,1];
I_vec=2e8*[1,1,1];

% generation of coordinates and connectivities
x1=3000*(1+cos(pi/4));
nodeCoordinates=[0 3000;
                 3000 3000;
                 x1 0;
                 x1+3000 0];
elementNodes=[1 2;
              2 3;
              3 4];
numberNodes=size(nodeCoordinates,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% Assembly stiffness matrix
K_Assembly=formStiffness2Dframe(GDof,size(elementNodes,1),elementNodes,nodeCoordinates,E_vec,I_vec,A_vec);

%force vector
F_col=nan(GDof,1);
F_col(4)=0;
F_col(5)=-10000;
F_col(6)=-5e6;
F_col(7)=0;
F_col(8)=-10000;
F_col(9)=5e6;

%displacement vector
prescribedDof=[1 2 3 10 11 12];
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_Assembly,D_col,F_col);

% % drawing undeformed and deformed meshes
% U=D_col;  
% drawingMesh(nodeCoordinates+500*[U(1:numberNodes)...
%     U(numberNodes+1:2*numberNodes)],elementNodes,'L2','k.-');
% drawingMesh(nodeCoordinates,elementNodes,'L2','k--');
