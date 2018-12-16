function [nodeCoordinates, elementNodes]=rectangularMesh(L_x,L_y,N_elements_x,N_elements_y)

N_elements=N_elements_x*N_elements_y;
N_nodes_x=N_elements_x+1;
N_nodes_y=N_elements_y+1;

x_vec=linspace(0,L_x,N_nodes_x);
y_vec=linspace(0,L_y,N_nodes_y);
[x_mat,y_mat]=meshgrid(x_vec,y_vec);


nodeCoordinates=[x_mat(:),y_mat(:)];
elementNodes=nan(N_elements,4);
for ii=1:N_elements_x
    iStart=(ii-1)*N_nodes_y+1;
    ind_sw=iStart+(1:N_elements_y).'-1;
    ind_nw=ind_sw+1;
    ind_ne=ind_nw+N_nodes_y;
    ind_se=ind_sw+N_nodes_y;
    elementNodes((ii-1)*N_elements_y+(1:N_elements_y),:)=[ind_se ind_ne ind_nw ind_sw];
end

