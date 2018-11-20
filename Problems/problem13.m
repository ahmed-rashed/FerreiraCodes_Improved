% antonio ferreira 2008
%Modified by Ahmed Rashed
%EOV in this example are not clear

close all

numberElements=8;
E_vec=210e6*ones(1,numberElements);
G_vec=84e6*ones(1,numberElements);

A_vec=0.02*ones(1,numberElements);
Iy_vec=10e-5*ones(1,numberElements);
Iz_vec=20e-5*ones(1,numberElements);
J_vec=5e-5*ones(1,numberElements);

% generation of coordinates and connectivities
nodesCoords=[   0 0 0; 
                0 0 4; 
                4 0 4; 
                4 0 0;
                0 5 0; 
                0 5 4; 
                4 5 4; 
                4 5 0];

elementNodes=[  1 5
                2 6
                3 7
                4 8
                5 6
                6 7
                7 8
                8 5]; 
numberNodes=size(nodesCoords,1);

% GDof: global number of degrees of freedom
GDof=6*numberNodes; 

% calculation of the system stiffness matrix
K_assembly=formStiffness3Dframe(GDof,numberElements,elementNodes,nodesCoords,E_vec,A_vec,Iz_vec,Iy_vec,G_vec,J_vec);

% boundary conditions and solution
prescribedDof=1:4*6;

%force vector
F_col=nan(GDof,1);
F_col(4*6+1:end)=0;
F_col(37)=-15;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_assembly,D_col,F_col);
