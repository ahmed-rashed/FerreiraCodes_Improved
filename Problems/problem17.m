% 2D problem: thin plate in tension
% antonio ferreira 2008
% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

clc
clearvars
close all

% materials
E=10e7;
nu=0.30;  

%Plane stress model
C=E/(1-nu^2)*[  1 nu 0
                nu 1 0
                0 0 (1-nu)/2];

% load
P = 1e6;

%Mesh generation
L_x=5;
L_y=1;
N_elements_X=20;
N_elements_Y=10;
N_elements=N_elements_X*N_elements_Y;
%Create your "rectangularMesh"
[nodeCoordinates, elementNodes]=rectangularMesh(L_x,L_y,N_elements_X,N_elements_Y);
% Replace the following command
%drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
N_nodes=size(nodeCoordinates,1);

% GDof: global number of degrees of freedom
GDof=2*N_nodes; 

% calculation of the system stiffness matrix
[K_Assembly,M_Assembly]=formStiffness2D(GDof,N_elements,elementNodes,nodeCoordinates,C,1,1);

% boundary conditions 
iNodeLeftEdge=find(nodeCoordinates(:,1)==0);
iNodeBottomEdge=find(nodeCoordinates(:,2)==0);
iNodeRightEdge=find(nodeCoordinates(:,1)==L_x);
iNodeRightTop=iNodeRightEdge(nodeCoordinates(iNodeRightEdge,2)==L_y);
iNodeRightbottom=iNodeRightEdge(nodeCoordinates(iNodeRightEdge,2)==0);

prescribedDof=[ 2*iNodeLeftEdge-1
                2*iNodeBottomEdge];

% force vector (distributed load applied at x=L_x)
F_col=zeros(GDof,1);
F_col(prescribedDof)=nan;

F_equiv=zeros(GDof,1);
F_equiv(2*iNodeRightEdge-1)=P*(L_y/N_elements_Y);
F_equiv(2*iNodeRightTop-1)=P*(L_y/N_elements_Y)/2;
F_equiv(2*iNodeRightbottom-1)=P*(L_y/N_elements_Y)/2;

%displacement vector
D_col=nan(GDof,1);
D_col(prescribedDof)=0;

% solution
[D_col,F_col]=solution(prescribedDof,K_Assembly,D_col,F_col,F_equiv);

%Normal Modes Analysis
N_modes=4;
[D_modeShape_cols,w_n_vec]=solutionModal(prescribedDof,D_col(prescribedDof),K_Assembly,M_Assembly,N_modes);

% Drawing
matrixShape=[N_elements_Y+1,N_elements_X+1];
scaleFactor=10;

D_x_mat=reshape(D_col(1:2:end),matrixShape);
D_y_mat=reshape(D_col(2:2:end),matrixShape);
x_mat=reshape(nodeCoordinates(:,1),matrixShape);
y_mat=reshape(nodeCoordinates(:,2),matrixShape);

x_deformed_mat=x_mat+scaleFactor*D_x_mat;
y_deformed_mat=y_mat+scaleFactor*D_y_mat;
D_magnitude_mat=sqrt(D_x_mat.^2+D_y_mat.^2);
h=pcolor(x_deformed_mat,y_deformed_mat,D_x_mat);
set(h,'EdgeColor',.5*[1,1,1])
hold on

mesh(x_mat,y_mat,0*x_mat,'FaceColor','none','EdgeColor','k')
view(2)
axis equal
axis off
title('Deformed shape')
colorbar

% stresses at nodes
stresses=stresses2D(N_elements,elementNodes,N_nodes,nodeCoordinates,D_col,D_x_mat(:),D_y_mat(:),C,scaleFactor);