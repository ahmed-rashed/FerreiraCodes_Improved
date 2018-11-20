function [AMatrix,BMatrix,DMatrix,SMatrix,qbarra]=...
    srinivasMaterial(thickness)

%%% SRINIVAS EXAMPLE

dd=zeros(2);d=zeros(3);
% multiplying factor for skins
rf=15;
% plate thickness
h=thickness;
% matrix [D]
% in-plane
dmat=[0.999781 0.231192 0 ;0.231192 0.524886 0;0 0 0.262931];
% shear
dm=[0.26681 0; 0 0.159914];
% nc: number of layers
nc=3; 
% layers angles
ttt=0;ttt1=0; th(1)=ttt;th(2)=ttt1;th(3)=ttt1;
% coordinates - z1 (upper) and - z2 (lower) for each layer
z1=[-2*h/5 2*h/5 h/2];
z2=[-h/2 -2*h/5 2*h/5;];
% thickness for each layer
thick(1:nc)=z1(1:nc)-z2(1:nc);
% coefe: shear correction factors (k1 and k2)
coefe(1:2)=0.0;gbarf(1:2)=0.0;rfact(1:2)=0.0;
sumla(1:2)=0.0;trlow(1:2)=0.0;upter(1:2)=0.0;
% middle axis position (bending)
      dsumm=0.0;
      for ilayr=1:nc
        dzeta=z1(ilayr)-z2(ilayr);         
        zheig=dsumm+dzeta/2.0;
                 
            dindx(1)=rf*dmat(1,1);dindx(2)=dmat(2,2);
            upter(1:2)=upter(1:2)+dindx(1:2)*zheig*dzeta;            
            trlow(1:2)=trlow(1:2)+dindx(1:2)*dzeta;
                                   
         dsumm = dsumm+dzeta;
      end
     
      zeta2(1:2)=-upter(1:2)./trlow(1:2);

% shear correction factors.

      for  ilayr=1:nc
            diff1=z1(ilayr)-z2(ilayr);
            d1=rf*dmat(1,1);
            d2=rf*dmat(2,2);
            d3=rf*dm(1,1);
            d4=rf*dm(2,2);
            if(ilayr==2)
                d1=dmat(1,1);
                d3=dm(1,1);
                d4=dm(2,2);
                d2=dmat(2,2);
            end
         index=10;
         for i=1:2
            zeta1(i)=zeta2(i);
            zeta2(i)=zeta1(i)+diff1;
            diff2(i)=zeta2(i)^2-zeta1(i)^2;
            diff3(i)=zeta2(i)^3-zeta1(i)^3;
            diff5(i)=zeta2(i)^5-zeta1(i)^5;

            dindx=[d1;d2];      
            gindx=[d3;d4];    
           
            gbarf(i)=gbarf(i)+gindx(i)*diff1/2.0;     
            rfact(i)=rfact(i)+dindx(i)*diff3(i)/3.0;

            term1   = sumla(i)*sumla(i)*diff1;       
            term2   = dindx(i)*(zeta1(i)^4)*diff1/4.0;
            term3   = dindx(i)*diff5(i)/20.0;         
            term4   =-dindx(i)*zeta1(i)*zeta1(i)*diff3(i)/6.0;
            term5   = sumla(i)*zeta1(i)*zeta1(i)*diff1;    
            term6   =-sumla(i)*diff3(i)/3.0;
            coefe(i)= coefe(i)+(term1+dindx(i)*...
                (term2+term3+term4+term5+term6))/gindx(i);    
            index   = index+1;
            sumla(i)= sumla(i)-dindx(i)*diff2(i)/2.0;
        end
      end
      
    coefe(1:2)=rfact(1:2).*rfact(1:2)./(2.0*gbarf(1:2).*coefe(1:2));
    kapa=coefe(1);
    
% constitutive matrice, membrane, bending and shear
a11=0;a22=0;a12=0;a33=0;
for i=1:nc
theta=th(i);
q11=rf*dmat(1,1);q12=rf*dmat(1,2);q22=rf*dmat(2,2);q33=rf*dmat(3,3);
cs=cos(theta);ss=sin(theta);ss11=rf*dm(1,1)*kapa;ss22=rf*dm(2,2)*kapa;
if i==2
    q11=dmat(1,1);q12=dmat(1,2);q22=dmat(2,2);q33=dmat(3,3);
    cs=cos(theta);ss=sin(theta);
    ss11=dm(1,1)*kapa;ss22=dm(2,2)*kapa;
end
dd(1,1)=dd(1,1)+(ss11*cos(theta)^2+ss22*sin(theta)^2)*(z1(i)-z2(i));
dd(2,2)=dd(2,2)+(ss11*sin(theta)^2+ss22*cos(theta)^2)*(z1(i)-z2(i));
d(1,1)=d(1,1)+(q11*cs^4+2*(q12+2*q33)*ss*ss*cs*cs+...
    q22*ss^4)*(z1(i)^3-z2(i)^3)/3;
d(2,2)=d(2,2)+(q11*ss^4+2*(q12+2*q33)*ss*ss*cs*cs+...
    q22*cs^4)*(z1(i)^3-z2(i)^3)/3;
d(1,2)=d(1,2)+((q11+q22-4*q33)*ss*ss*cs*cs+...
    q12*(ss^4+cs^4))*(z1(i)^3-z2(i)^3)/3;
d(3,3)=d(3,3)+((q11+q22-2*q12-2*q33)*ss*ss*cs*cs+...
    q33*(ss^4+cs^4))*(z1(i)^3-z2(i)^3)/3;
a11=a11+q11*thick(i);
a22=a22+q22*thick(i);
a33=a22+q33*thick(i);
a12=a12+q12*thick(i);

qbarra(1,1,i)=q11;
qbarra(1,2,i)=q12;
qbarra(2,2,i)=q22;
qbarra(3,3,i)=q33;
qbarra(4,4,i)=ss11;
qbarra(5,5,i)=ss22;

end %nc

A44=dd(2,2);
A55=dd(1,1);
D11=d(1,1);
D12=d(1,2);
D22=d(2,2);
D66=d(3,3);
A11=a11;
A12=a12;
A66=a33;
A22=a22;

AMatrix=[A11,A12,0;A12,A22,0;0,0,A66];
%srinivas case (symmetric)
BMatrix=zeros(3);
%BMatrix=[B11,B12,0;B12,B22,0;0,0,B66]
DMatrix=[D11,D12,0;D12,D22,0;0,0,D66];
SMatrix=[A44,0;0,A55];