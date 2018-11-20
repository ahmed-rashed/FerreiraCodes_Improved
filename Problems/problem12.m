% antonio ferreira 2008
%Modified by Ahmed Rashed
%EOV in this example are not clear

close all
clearvars

E_vec=210e6*[1,1,1];
G_vec=84e6*[1,1,1];

A_vec=0.02*[1,1,1];
Iy_vec=10e-5*[1,1,1];
Iz_vec=20e-5*[1,1,1];
J_vec=5e-5*[1,1,1];

% generation of coordinates and connectivities
nodesCoords=[   0 0 0;
                3 0 0;
                0 0 -3;
                0 -4 0];

elementNodes=[  1 2;
                1 3;
                1 4];
numberNodes=size(nodesCoords,1);

% GDof: global number of degrees of freedom
GDof=6*numberNodes;

% Assembly stiffness matrix
K_assembly=formStiffness3Dframe(GDof,size(elementNodes,1),elementNodes,nodesCoords,E_vec,A_vec,Iz_vec,Iy_vec,G_vec,J_vec);

% boundary conditions
prescribedDof=7:24;

%force vector
F_col=nan(GDof,1);
F_col([2,4:6])=0;
F_col(1)=-10;
F_col(3)=20;    %This is diferent from Figure 8.1

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);
