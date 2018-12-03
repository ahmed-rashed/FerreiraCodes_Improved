function  k_local=k_local_3D_Frame(L,E,A,Iz,Iy,G,J)
a=[E*A/L,0,0
    0,12*E*Iz/L^3,0
    0,0,12*E*Iy/L^3];

b=[ 0,0,0
    0,0,6*E*Iz/L^2
    0,-6*E*Iy/L^2,0];

bb=[ 0,0,0
    0,0,6*E*Iy/L^2
    0,-6*E*Iz/L^2,0];

c=[G*J/L,0,0
    0,4*E*Iy/L,0
    0,0,4*E*Iz/L];

d=[-G*J/L,0,0
    0,2*E*Iy/L,0
    0,0,2*E*Iz/L];

k_local=[a,b,-a,b;
       b.',c,bb,d;
        -a.',bb.',a,-b;
        b.',d.',-b.',c];