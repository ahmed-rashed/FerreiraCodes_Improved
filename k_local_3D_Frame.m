function  k_local=k_local_3D_Frame(nodeCoordinates,iNodes,E,A,Iz,Iy,G,J)

L=norm(nodeCoordinates(iNodes(2),:)-nodeCoordinates(iNodes(1),:));

k1 = E*A/L;
k2 = 12*E*Iz/(L^3);
k3 = 6*E*Iz/(L^2);
k4 = 4*E*Iz/L;
k5 = 2*E*Iz/L;
k6 = 12*E*Iy/(L^3);
k7 = 6*E*Iy/(L^2);
k8 = 4*E*Iy/L;
k9 = 2*E*Iy/L;
k10 = G*J/L;

a=[k1 0 0; 0 k2 0; 0 0 k6];
b=[ 0 0 0;0 0 k3; 0 -k7 0];
c=[k10 0 0;0 k8 0; 0 0 k4];
d=[-k10 0 0;0 k9 0;0 0 k5];

k_local=[a b -a b;
        b.' c b d;
        -a.' b.' a -b;
        b.' d.' -b.' c];