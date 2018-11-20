clear all;
% 
E0  = 10e6;  % m—dulo E
nu0 = 0.333;  % Poisson
%
L  = 1;     % comprimento 
thickness=0.001;
rho=1;
I=thickness^3/12;% momento de inercia
A=1*thickness;
% 
numx     = 40;  
%
P = 1; % pressao uniforme
% matriz constitutiva
G=E0/2/(1+nu0);
C=[   I*E0     0; 0    5/6*thickness*G];
%
% malha
node=linspace(0,L,numx+1);xx=node;x=node;
%
for i=1:size(node,2)-1
    element(i,1)=i; 
    element(i,2)=i+1
end
%  
numnode=size(node,2);    % num. nos
numelem=size(element,1); % num. elementos
%
K=zeros(2*numnode,2*numnode);
Kg=zeros(2*numnode,2*numnode);
% rigidez
W=zeros(2); Q=zeros(2);

Q(1) = 0.577350269189626;
Q(2) =-0.577350269189626;

W(1) = 1.;  W(2) = 1.; 

for e=1:numelem
    indice=element(e,:);  
    indiceB=[ indice indice+numnode]; 
    indiceR=indice+numnode;   
    nn=length(indice);    
    length_element=xx(indice(2))-xx(indice(1));
    detJ0=length_element/2;invJ0=1/detJ0;
  for q=1:size(W,1) ; 
     pt=Q(q,:); wt=W(q);  
     pt=pt(1);
     N=([1-pt,1+pt]/2)';
     dNdxi=[-1;1]/2;
     dNdx=dNdxi*invJ0;
% B     
B=zeros(2,2*nn);  B(1,nn+1:2*nn)  = dNdx(:)'; 
% K
K(indiceB,indiceB)=K(indiceB,indiceB)+B'*B*W(q)*detJ0*C(1,1);
Kg(indice,indice)=Kg(indice,indice)+dNdx*dNdx'*W(q)*detJ0*P;
  end  
end   
%
W=zeros(1); Q=zeros(1);
Q(1) = 0.;
W(1) = 2.;  

for e=1:numelem  
    indice=element(e,:); 
    indiceB=[ indice indice+numnode];
    nn=length(indice);  
    length_element=xx(indice(2))-xx(indice(1));
    detJ0=length_element/2;invJ0=1/detJ0;
  for q=1:size(W,1) ;     
     pt=Q(q,:); wt=W(q);  
     pt=pt(1);
     N=([1-pt,1+pt]/2)';
     dNdxi=[-1;1]/2;
     dNdx=dNdxi*invJ0;
% B     
    B=zeros(2,2*nn);    
    B(2,1:nn)       = dNdx(:)';  
    B(2,nn+1:2*nn)  = N;
% K
    K(indiceB,indiceB)=K(indiceB,indiceB)+B'*B*W(q)*detJ0*C(2,2); 
  end  
end    
% BC 
%
% CC
fixedNodeW =find(xx==min(node(:)) | xx==max(node(:)))';
fixedNodeTX=fixedNodeW; 

% CF
%fixedNodeW =find(xx==min(node(:)))';
%fixedNodeTX=fixedNodeW;         

% SS
fixedNodeW =find(xx==min(node(:)) | xx==max(node(:)))';
fixedNodeTX=[];         
% 
dofs=[fixedNodeW;fixedNodeTX+numnode]                           
activeDof=setdiff([1:2*numnode]',[dofs]);

[V,D]=eig(K(activeDof,activeDof),Kg(activeDof,activeDof));
D=diag(D);[D,ii] = sort(D);  V = V(:,ii);

kapa=5/6;
PcrSS=pi*pi*E0*I/L^2*(1/(1+pi*pi*E0*I/(L*L*kapa*G*A)))
PcrCC=pi*pi*E0*I/(L/2)^2*(1/(1+pi*pi*E0*I/(L*L/4*kapa*G*A)))
    
    modeNumber=4;

    V1=zeros(2*numnode,modeNumber);
    V1(activeDof,1:modeNumber)=V(:,1:modeNumber);
        
    clf
    for j=1:modeNumber
    u=[V1(1:numnode,j)];
    xx = 0:.01:1; subplot(modeNumber+1,1,j);
    plot(x',u,'.','markersize',12), grid on
    uu = polyval(polyfit(x',V1(1:numnode,j),numnode),xx);
    line(xx,uu)%, axis off
    end
    

   

