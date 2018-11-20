%................................................................

% MATLAB codes for Finite Element Analysis
% problem20.m
% laminated plate: Srinivas problem:
% S. Srinivas, A refined analysis of composite laminates, 
% J. Sound and Vibration, 30 (1973),495--507.

% antonio ferreira 2008


colordef white
% propriedades material
E0  = 10920;  % m—dulo Young
nu0 = 0.30;  % coef. Poisson
rho=1;
% geometria da placa de lado L
L  = 1;     % lado
thickness=0.1;
I=thickness^3/12;
%carga
P=-1;
numy     = 20;    % # elementos em xx 
numx     = 20;    % # elementos em xx 

%%% SRINIVAS EXAMPLE

dd=zeros(2);d=zeros(3);
% factor de multiplicacao de propriedades materiais
rf=15;
% espssura da placa
h=thickness;
% matrix [D]
% plano
dmat=[0.999781 0.231192 0 ;0.231192 0.524886 0;0 0 0.262931];
% corte
dm=[0.26681 0; 0 0.159914];
% nc: numero de camadas
nc=3; 
% angulos das camadas
ttt=0;ttt1=0; 
th(1)=ttt;th(2)=ttt1;th(3)=ttt1;
% coordenadas z1 (superior) e z2 (inferior) de cada camada
z1=[-2*h/5 2*h/5 h/2];
z2=[-h/2 -2*h/5 2*h/5;];
% espessura de cada camada
thick(1:nc)=z1(1:nc)-z2(1:nc);
% coefe: factores de correccao ao corte (k1 e k2)
coefe(1:2)=0.0;gbarf(1:2)=0.0;rfact(1:2)=0.0;
sumla(1:2)=0.0;trlow(1:2)=0.0;upter(1:2)=0.0;
% posicao do eixo neutro de flexao
      dsumm=0.0;
      for ilayr=1:nc
        dzeta=z1(ilayr)-z2(ilayr);         
        zheig=dsumm+dzeta/2.0;
         for i=1:2
            dindx(i)=rf*dmat(i,i);
            if(ilayr==2)
                 dindx(i)=dmat(i,i);
            end
            upter(i)=upter(i)+dindx(i)*zheig*dzeta;            
            trlow(i)=trlow(i)+dindx(i)*dzeta;
        end
         dsumm = dsumm+dzeta;
      end
     
      zeta2(1:2)=-upter(1:2)./trlow(1:2);

% calculo dos factores de correccao ao corte.

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
            if(i==1)
                dindx(i)=d1;      
                gindx(i)=d3;
            else
                dindx(i)=d2; 
                gindx(i)=d4;
            end                
           
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
% calculo das matrizes constitutivas de membrana, flexao e corte
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

AMatrix=[A11,A12,0;A12,A22,0;0,0,A66]
%srinivas case (symmetric)
BMatrix=zeros(3);
%BMatrix=[B11,B12,0;B12,B22,0;0,0,B66]
DMatrix=[D11,D12,0;D12,D22,0;0,0,D66]
SMatrix=[A44,0;0,A55]

% geracao de malha
[node, element] = MalhaRectangular(L, L, numx, numy);
%
xx=node(:,1);   yy=node(:,2);
numnode=size(node,1);    % numero nos
numelem=size(element,1); % numero elementos
% inicializacao de matrizes e vectores
K=zeros(5*numnode,5*numnode); % matriz de rigidez
M=zeros(5*numnode,5*numnode); % matriz de massa
f=zeros(5*numnode,1); % vector de carga
U=zeros(5*numnode,1); % vector de deslocamentos

% calculo de MATRIZ DE RIGIDEZ

% quadratura de 2 pontos
Q=[ -0.577350269189626 -0.577350269189626;
    0.577350269189626 -0.577350269189626;
    0.577350269189626 0.577350269189626;
    -0.577350269189626 0.577350269189626];
W=[ 1;1;1;1]; 
% rigidez de flexao 
for e=1:numelem                          
  % indice: nos do elemento (e)
  indice=element(e,:);    
  % indiceB: vector com graus de liberdade do elemento (e)
  indiceB=[ indice indice+numnode indice+2*numnode ...
            indice+3*numnode indice+4*numnode ];    
  nn=length(indice);
  for q=1:size(W,1)  % CICLO QUADRATURA
    pt=Q(q,:);       % PONTO de Gauss
    wt=W(q);         % PESO de Gauss
    xi=pt(1);        % coordenadas naturais, xi, eta
    eta=pt(2);
    % FUNCOES DE FORMA E DERIVADAS COORDENADAS NATURAIS
    N=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
            (1+xi)*(1+eta);(1-xi)*(1+eta)];
    dNdxi=1/4*[-(1-eta),     -(1-xi);1-eta,    -(1+xi);
		         1+eta,      1+xi;-(1+eta),   1-xi]; 
    % MATRIZ JACOBIANA 
    J0=node(indice,:)'*dNdxi;        
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    
% MATRIZ B   : flexao .....
    B_b=zeros(3,5*nn);
    B_b(1,nn+1:2*nn)        = dNdx(:,1)';  
    B_b(2,2*nn+1:3*nn)      = dNdx(:,2)';
    B_b(3,nn+1:2*nn)        = dNdx(:,2)';  
    B_b(3,2*nn+1:3*nn)      = dNdx(:,1)';
% MATRIZ B   : membrana .....
    B_m=zeros(3,5*nn);
    B_m(1,3*nn+1:4*nn)      = dNdx(:,1)';  
    B_m(2,4*nn+1:5*nn)      = dNdx(:,2)';
    B_m(3,3*nn+1:4*nn)      = dNdx(:,2)';  
    B_m(3,4*nn+1:5*nn)      = dNdx(:,1)';

% RIGIDEZ NO PONTO DE QUADRATURA
% note-se que aqui ja se faz o espalhamento da matriz

% ... flexao-flexao
    K(indiceB,indiceB)=K(indiceB,indiceB)+...
                       B_b'*DMatrix*B_b*W(q)*det(J0);
% ... membrana-membrana                  
    K(indiceB,indiceB)=K(indiceB,indiceB)+...
                       B_m'*AMatrix*B_m*W(q)*det(J0);
% ... membrana-flexao                  
    K(indiceB,indiceB)=K(indiceB,indiceB)+...
                       B_m'*BMatrix*B_b*W(q)*det(J0);
% ... flexao-membrana                  
    K(indiceB,indiceB)=K(indiceB,indiceB)+...
                       B_b'*BMatrix*B_m*W(q)*det(J0);

% VECTOR DE CARGA NO PONTO DE QUADRATURA
    f(indice)=f(indice)+N*P*det(J0)*wt;    
  end  
end    

%  MATRIZ DE RIGIDEZ (CORTE)
% quadratura de 1 ponto
 Q=[0 0];
 W=[4];

for e=1:numelem                         
  % indice: nos do elemento (e)
  indice=element(e,:);    
  % indiceB: vector com graus de liberdade do elemento (e)
  indiceB=[ indice indice+numnode indice+2*numnode ...
            indice+3*numnode indice+4*numnode];    
  nn=length(indice);
  for q=1:size(W,1)  % CICLO QUADRATURA
    pt=Q(q,:);       % PONTO de Gauss
    wt=W(q);         % PESO de Gauss
    xi=pt(1);        % coordenadas naturais, xi, eta
    eta=pt(2);
    % FUNCOES DE FORMA E DERIVADAS COORDENADAS NATURAIS
    N=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
            (1+xi)*(1+eta);(1-xi)*(1+eta)];
    dNdxi=1/4*[-(1-eta),     -(1-xi);1-eta,    -(1+xi);
		         1+eta,      1+xi;-(1+eta),   1-xi]; 
    % MATRIZ JACOBIANA 
    J0=node(indice,:)'*dNdxi;        
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    
% MATRIZ B   : corte ..... 
    B_s=zeros(2,5*nn);
    B_s(1,1:nn)       = dNdx(:,1)';  
    B_s(2,1:nn)       = dNdx(:,2)';
    B_s(1,nn+1:2*nn)  = N;           
    B_s(2,2*nn+1:3*nn)= N;

% ... corte-corte
    K(indiceB,indiceB)=K(indiceB,indiceB)+...
        B_s'*SMatrix*B_s*W(q)*det(J0);
  end  
end    

% imposicao de condicoes fronteira essenciais
% nos com w restringido (SSSS ou CCCC)
fixedNodeW =find(yy==max(node(:,2))|xx==min(node(:,1))...
                |xx==max(node(:,1))|yy==min(node(:,2)));
% caso simplesmente apoiado
fixedNodeTX =find(yy==max(node(:,2))|yy==min(node(:,2)));
fixedNodeTY =find(xx==max(node(:,1))| xx==min(node(:,1)));
fixedNodeU =find(xx==min(node(:,1)));
fixedNodeV =find(yy==min(node(:,2)));

% caso encastrado: nota , descomentar as proximas linhas e
% comentar as anteriores, conforme o caso que se quiser testar
% fixedNodeTX =fixedNodeW;
% fixedNodeTY =fixedNodeTX;
% fixedNodeU =fixedNodeTX;
% fixedNodeV =fixedNodeTX;

% vector com todas as condicoes fronteira
dofs=[fixedNodeW;fixedNodeTX+numnode;...
      fixedNodeTY+2*numnode;...
      fixedNodeU+3*numnode;fixedNodeV+4*numnode];
% graus de liberdade activos
activeDof=setdiff([1:5*numnode]',[dofs]);
% solucao do sistema de equacoes
U=K([activeDof],[activeDof])\f([activeDof]);
U1=zeros(5*numnode,1);
U1(activeDof)=U;
% U: vector solucao final
U=U1;
% desenho de deformada
ws=1:numnode;
% normalizacao para solucao de SRINIVAS
min(U1(ws))*0.999781/h
figure (1)
plot3(xx,yy,U(ws),'.')

% tensoes normais
% em cada camada
  stress_camada1=zeros(numelem,4,3);
  stress_camada2=zeros(numelem,4,3);
  stress_camada3=zeros(numelem,4,3);
% tensoes normais para cada camada
  stress_camada1=zeros(numelem,4,3);
  stress_camada2=zeros(numelem,4,3);
  stress_camada3=zeros(numelem,4,3);
% quadratura de 2 pontos
Q=[ -0.577350269189626 -0.577350269189626;
    0.577350269189626 -0.577350269189626;
    0.577350269189626 0.577350269189626;
    -0.577350269189626 0.577350269189626];
W=[ 1;1;1;1]; 

for e=1:numelem                         
  % indice: nos do elemento (e)
  indice=element(e,:);    
  % indiceB: vector com graus de liberdade do elemento (e)
  indiceB=[ indice indice+numnode indice+2*numnode ...
            indice+3*numnode indice+4*numnode];    
  nn=length(indice);
  for q=1:size(W,1)  % CICLO QUADRATURA
    pt=Q(q,:);       % PONTO de Gauss
    wt=W(q);         % PESO de Gauss
    xi=pt(1);        % coordenadas naturais, xi, eta
    eta=pt(2);
    % FUNCOES DE FORMA E DERIVADAS COORDENADAS NATURAIS
    N=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
            (1+xi)*(1+eta);(1-xi)*(1+eta)];
    dNdxi=1/4*[-(1-eta),     -(1-xi);1-eta,    -(1+xi);
		         1+eta,      1+xi;-(1+eta),   1-xi]; 
    % MATRIZ JACOBIANA 
    J0=node(indice,:)'*dNdxi;        
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    
% MATRIZ B   : flexao .....    
    B_b=zeros(3,5*nn);
    B_b(1,nn+1:2*nn)        = dNdx(:,1)';  
    B_b(2,2*nn+1:3*nn)      = dNdx(:,2)';
    B_b(3,nn+1:2*nn)        = dNdx(:,2)';  
    B_b(3,2*nn+1:3*nn)      = dNdx(:,1)';
% MATRIZ B   : membrana .....        
    B_m=zeros(3,5*nn);
    B_m(1,3*nn+1:4*nn)      = dNdx(:,1)';  
    B_m(2,4*nn+1:5*nn)      = dNdx(:,2)';
    B_m(3,3*nn+1:4*nn)      = dNdx(:,2)';  
    B_m(3,4*nn+1:5*nn)      = dNdx(:,1)';   
    
    stress_camada1(e,q,:)=...
        2*h/5*qbarra(1:3,1:3,2)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,2)*B_m*U(indiceB);
    stress_camada2(e,q,:)=...
        2*h/5*qbarra(1:3,1:3,3)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,3)*B_m*U(indiceB);
    stress_camada3(e,q,:)=...
        h/2*qbarra(1:3,1:3,3)*B_b*U(indiceB)+...
        qbarra(1:3,1:3,3)*B_m*U(indiceB);
  end                               
end   

% tensoes de corte em cada camada
% directamente pelas equacoes constitutivas
% do material de cada camada

    shear_camada1=zeros(numelem,1,2);
    shear_camada2=zeros(numelem,1,2);
    shear_camada3=zeros(numelem,1,2);

% quadratura de 1 ponto
 Q=[0 0];
 W=[4];


for e=1:numelem                         
  
  % indice: nos do elemento (e)
  indice=element(e,:);    
  % indiceB: vector com graus de liberdade do elemento (e)
  indiceB=[ indice indice+numnode indice+2*numnode ...
            indice+3*numnode indice+4*numnode];    
  nn=length(indice);
  for q=1:size(W,1)  % CICLO QUADRATURA
    pt=Q(q,:);       % PONTO de Gauss
    wt=W(q);         % PESO de Gauss
    xi=pt(1);        % coordenadas naturais, xi, eta
    eta=pt(2);
    % FUNCOES DE FORMA E DERIVADAS COORDENADAS NATURAIS
    N=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
            (1+xi)*(1+eta);(1-xi)*(1+eta)];
    dNdxi=1/4*[-(1-eta),     -(1-xi);1-eta,    -(1+xi);
		         1+eta,      1+xi;-(1+eta),   1-xi]; 
    % MATRIZ JACOBIANA 
    J0=node(indice,:)'*dNdxi;        
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    
% CORTE   
    B_s=zeros(2,5*nn);
    B_s(1,1:nn)       = dNdx(:,1)';  
    B_s(2,1:nn)       = dNdx(:,2)';
    B_s(1,nn+1:2*nn)  = N;           
    B_s(2,2*nn+1:3*nn)= N;
            
    shear_camada1(e,q,:)=qbarra(4:5,4:5,1)*B_s*U(indiceB);
    shear_camada2(e,q,:)=qbarra(4:5,4:5,2)*B_s*U(indiceB);

  end                               
end  
% apresentacao das tensoes, ja normalizadas
% para o problema de Srinivas
format
[ abs(min(U1(ws))*0.999781/thickness), ...
  abs(min(stress_camada3(:,3,1))),...
  abs(min(stress_camada2(:,3,1))), ...
  abs(min(stress_camada1(:,3,1))),...
  max(shear_camada2(:,:,1)),...
  max(shear_camada1(:,:,1))]
