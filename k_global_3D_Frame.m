function  k_global=k_global_3D_Frame(nodeCoordinates,iNodes,k_local)

L=norm(nodeCoordinates(iNodes(2),:)-nodeCoordinates(iNodes(1),:));

x1=nodeCoordinates(iNodes(1),1);
y1=nodeCoordinates(iNodes(1),2);
z1=nodeCoordinates(iNodes(1),3);
x2=nodeCoordinates(iNodes(2),1);
y2=nodeCoordinates(iNodes(2),2);
z2=nodeCoordinates(iNodes(2),3);

if x1 == x2 && y1 == y2
   if z2 > z1
      T = [0 0 1;
           0 1 0;
          -1 0 0];
   else
      T = [0 0 -1;
           0 1 0;
           1 0 0];
   end
else
    CXx = (x2-x1)/L;
    CYx = (y2-y1)/L;
    CZx = (z2-z1)/L;
    D = sqrt(CXx^2 + CYx^2);
    CXy = -CYx/D;
    CYy = CXx/D;
    CZy = 0;
    CXz = -CXx*CZx/D;
    CYz = -CYx*CZx/D;
    CZz = D;
    T = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
end

TT=zeros(12,12);
for ii=1:3:10
    TT(ii:ii+2,ii:ii+2)=T;
end

k_global=TT.'*k_local*TT;
