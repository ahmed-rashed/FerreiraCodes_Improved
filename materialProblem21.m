function [CMembranaMembrana,CMembranaFlexao0,...
          CFlexao0Flexao0,CCorte0Corte0]=...
          materialProblem21(h)

% liew
e2=1;e1=40*e2;g23=0.5*e2;g13=0.6*e2;g12=g13;miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;
% e2=1;e1=1*e2;g23=e2/2.5;g13=g23;g12=g13;miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;

%alfas=[0];% 1 camada
%alfas=[0,pi/2];% 2 camadas
%alfas=[-pi/4,pi/4];% 2 camadas
%alfas=[0,pi/2,0,pi/2,0,pi/2,0,pi/2]; % 8 camadas
%alfas=[0,pi/2,pi/2,0];% 4 camadas
alfas=[0,pi/2,0];% 3 camadas

%z(1)=-(h/2);z(2)=(h/2);% 1 camadas equidistantes.
%z(1)=-(h/2);z(2)=0;z(3)=(h/2);% 2 camadas equidistantes.
%z(1)=-h/2;z(2)=-3*h/8;z(3)=-h/4;z(4)=-(h/8);z(5)=0;z(6)=-z(4);z(7)=-z(3);z(8)=-z(2);z(9)=-z(1);% 8 camadas equidistantes.
%z(1)=-(h/2);z(2)=-(h/4);z(3)=0;z(4)=-z(2);z(5)=-z(1);%quatro camadas igualmente espaçadas
z(1)=-(h/2);z(2)=-(h/2)+h/3;z(3)=-z(2);z(4)=-z(1);%tres camadas igualmente espaçadas

%sem alteracao de angulo 0º 
qbarra(1,1,1)=e1/factor;
qbarra(1,2,1)=miu21*e1/factor;
qbarra(2,1,1)=miu12*e2/factor;
qbarra(2,2,1)=e2/factor;
qbarra(3,3,1)=g12;
qbarra(4,4,1)=kapa*g23;
qbarra(5,5,1)=kapa*g13;

T=[cos(phi)^2,sin(phi)^2,-sin(2*phi),0,0;...
   sin(phi)^2,cos(phi)^2,sin(2*phi),0,0;...
   sin(phi)*cos(phi),-sin(phi)*cos(phi),cos(phi)^2-sin(phi)^2,0,0;...
   0,0,0,cos(phi),sin(phi);...
   0,0,0,-sin(phi),cos(phi)];

qBarra=T*qbarra*T.';


for s=1:size(alfas,2)
    for i=1:5
        for j=1:5
           QQbarra(i,j,s)=subs(qBarra(i,j,1),phi,alfas(s));
       end
   end
   Qbarra=double(QQbarra);
end
Q=Qbarra; 

%______________________________________________
Astiff(5,5)=0;Bstiff(5,5)=0;Estiff(5,5)=0;Fstiff(5,5)=0;Gstiff(5,5)=0;Hstiff(5,5)=0;
Istiff(5,5)=0;Kstiff(5,5)=0;Lstiff(5,5)=0;
for k=1:size(alfas,2)
        for i=1:3
        for j=1:3
        Astiff(i,j)=Astiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
        Bstiff(i,j)=Bstiff(i,j)+Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
        Estiff(i,j)=Estiff(i,j)+Q(i,j,k)*(z(k+1)^4-z(k)^4)/4;
        Fstiff(i,j)=Fstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
        Gstiff(i,j)=Gstiff(i,j)+Q(i,j,k)*(z(k+1)^5-z(k)^5)/5; 
        Hstiff(i,j)=Hstiff(i,j)+Q(i,j,k)*(z(k+1)^7-z(k)^7)/7; 
        end
        end

            for i=4:5
            for j=4:5
            Istiff(i,j)=Istiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
            Kstiff(i,j)=Kstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
            Lstiff(i,j)=Lstiff(i,j)+Q(i,j,k)*(z(k+1)^5-z(k)^5)/5;
            end
            end
end
pi=double(pi);
%
CMembranaMembrana=Astiff(1:3,1:3);
CMembranaFlexao0=Bstiff(1:3,1:3);
CMembranaFlexao1=Estiff(1:3,1:3);
CFlexao0Flexao0=Fstiff(1:3,1:3);
CFlexao0Flexao1=Gstiff(1:3,1:3);
CFlexao1Flexao1=Hstiff(1:3,1:3);
CCorte0Corte0=Istiff(4:5,4:5);
CCorte0Corte1=Kstiff(4:5,4:5);
CCorte1Corte1=Lstiff(4:5,4:5);