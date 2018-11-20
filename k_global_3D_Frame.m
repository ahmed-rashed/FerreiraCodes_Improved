function  k_global=k_global_3D_Frame(nodeCoords,iNodes,k_local)

x1=nodeCoords(iNodes(1),1);
y1=nodeCoords(iNodes(1),2);
z1=nodeCoords(iNodes(1),3);
x2=nodeCoords(iNodes(2),1);
y2=nodeCoords(iNodes(2),2);
z2=nodeCoords(iNodes(2),3);

x_local=nodeCoords(iNodes(2),:)-nodeCoords(iNodes(1),:);
L=norm(x_local);
x_local_hat=x_local/L;  %Unit vector


D_x=x2-x1;
D_y=y2-y1;
D_z=z2-z1;

if x1 == x2 && y1 == y2
   if z2 > z1
      T_node = [0 0 1;
           0 1 0;
          -1 0 0];  %This assumes EOV is the global y axis
   else
      T_node = [0 0 -1;
           0 1 0;
           1 0 0];  %This assumes EOV is the global y axis
   end
else
    lx = D_x/L;
    mx = D_y/L;
    nx = D_z/L;
    
    D = sqrt(lx^2 + mx^2);
    ly = -mx/D;
    my = lx/D;
    ny = 0; %This assumes EOV is normal to the projection of x_local on the global xy plane
    
    lz = -lx*nx/D;
    mz = -mx*nx/D;
    nz = D;
    
    T_node = [x_local_hat
        ly my ny
        lz mz nz];
end

T_element=zeros(12,12);
for ii=1:3:10
    T_element(ii:ii+2,ii:ii+2)=T_node;
end

k_global=T_element.'*k_local*T_element;
