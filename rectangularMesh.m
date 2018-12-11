function [nodeCoordinates, elementNodes]=rectangularMesh(L_x,L_y,N_elements_X,N_elements_Y)

x_vec=linspace(0,L_x,N_elements_X);
y_vec=linspace(0,L_y,N_elements_Y);
[x_mat,y_mat]=meshgrid(x_vec,y_vec);

N_nodes=(N_elements_X+1)*(N_elements_Y+1);
N_elements=N_elements_X*N_elements_Y;

nodeCoordinates=[x_mat(:),y_mat(:)];
elementNodes=nan(N_elements,4);
=[]
