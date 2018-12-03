% 2D problem: thin plate in tension
% antonio ferreira 2008
%Modified by Ahmed Rashed

clearvars;
close all

% materials
E=10e7;
nu=0.30;  

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
fixedNodeX=find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY=find(nodeCoordinates(:,2)==0);  % fixed in YY
prescribedDof=[ fixedNodeX
                fixedNodeY+N_nodes];

% force vector (distributed load applied at xx=L_x)
force=zeros(GDof,1);
rightBord=find(nodeCoordinates(:,1)==L_x);
force(rightBord)   =P*L_y/N_elements_Y;
force(rightBord(1))=force(rightBord(1))/2;
force(rightBord(end))=force(rightBord(end))/2;

% solution
D_col=solution(GDof,prescribedDof,K_Assembly,force);

%Normal Modes Analysis
activeDof=setdiff(1:GDof,prescribedDof);

[modeShapes,lambda]=eig(K_Assembly(activeDof,activeDof),M_Assembly(activeDof,activeDof)); 
w_n=sqrt(lambda);

% sort out eigenvalues
[w_n,ii]=sort(w_n);
modeShapes=modeShapes(:,ii);

% % drawing eigenmodes
% N_modes=4;

% % D_col
% disp('Displacements')
% f=[1:GDof; D_col.'];
% fprintf('node U\n')
% fprintf('%3d %12.8f\n',f)
% UX=D_col(1:N_nodes);
% UY=D_col(N_nodes+1:GDof);
% scaleFactor=10;
% 
% % deformed shape
% figure
% drawingField(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4',UX);%U XX
% hold on
% drawingMesh(nodeCoordinates+scaleFactor*[UX UY],elementNodes,'Q4','k-');
% drawingMesh(nodeCoordinates,elementNodes,'Q4','k--');
% colorbar
% title('U XX  (on deformed shape)')
% axis off

% stresses at nodes
stresses=stresses2D(N_elements,elementNodes,N_nodes,nodeCoordinates,D_col,UX,UY,C,scaleFactor)