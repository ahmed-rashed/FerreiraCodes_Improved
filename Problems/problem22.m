%................................................................

% MATLAB codes for Finite Element Analysis
% problem22.m
% this function performs buckling analysis
% of Laminated plates using 5 degrees of freedom per node
% antonio ferreira 2008

function problem22

clearvars
colordef white

% L: side lenght
L  = 1;     

thickness=0.01;h=thickness;
I=thickness^3/12;
% kapa: shear correction factor
kapa=5/6;
%kapa=pi*pi/12;

% material properties
% symbolic computation 
syms phi

% liew material
e2=1;e1=40*e2;g23=0.5*e2;g13=0.6*e2;g12=g13;
miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;

% reddy material
e2=1;e1=25*e2;g23=0.2*e2;g13=0.5*e2;g12=g13;
miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;


% angles for laminate
alfas=[0,pi/2,pi/2,0];% 4 layers
% upper and lower coordinates
z(1)=-(h/2);z(2)=-(h/4);z(3)=0;z(4)=-z(2);z(5)=-z(1);
clear alfas
clear z
% angles for laminate
alfas=[0,pi/2,0];% 3 layers
% upper and lower coordinates
z(1)=-(h/2);z(2)=-(h/2)+h/3;z(3)=-z(2);z(4)=-z(1);


% [Q] in 0º orientation 
qbarra(1,1,1)=e1/factor;
qbarra(1,2,1)=miu21*e1/factor;
qbarra(2,1,1)=miu12*e2/factor;
qbarra(2,2,1)=e2/factor;
qbarra(3,3,1)=g12;
qbarra(4,4,1)=kapa*g23;
qbarra(5,5,1)=kapa*g13;

% transformation matrix
T=[cos(phi)^2,sin(phi)^2,-sin(2*phi),0,0;...
   sin(phi)^2,cos(phi)^2,sin(2*phi),0,0;...
   sin(phi)*cos(phi),-sin(phi)*cos(phi),cos(phi)^2-sin(phi)^2,0,0;...
   0,0,0,cos(phi),sin(phi);...
   0,0,0,-sin(phi),cos(phi)];

% [Q] in structural axes
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
Astiff(5,5)=0;Bstiff(5,5)=0;Fstiff(5,5)=0;Istiff(5,5)=0;
for k=1:size(alfas,2)
        for i=1:3
        for j=1:3
        Astiff(i,j)=Astiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
        Bstiff(i,j)=Bstiff(i,j)+Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
        Fstiff(i,j)=Fstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
        end
        end

            for i=4:5
            for j=4:5
            Istiff(i,j)=Istiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
            end
            end
end
%pi=double(pi); % come back to numeric computation

% constitutive matrices
CMembranaMembrana=Astiff(1:3,1:3);
CMembranaFlexao0=Bstiff(1:3,1:3);
CFlexao0Flexao0=Fstiff(1:3,1:3);
CCorte0Corte0=Istiff(4:5,4:5);

% initial stress matrix
sigmaX=1/thickness;
sigmaXY=0;
sigmaY=sigmaX;
sigmaMatrix=[ sigmaX sigmaXY; sigmaXY sigmaY];
                         
% mesh generation ...
% numberElementsX: number of elements in x
% numberElementsY: number of elements in y
numberElementsX=20;
numberElementsY=20;
% number of elements
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L, L, numberElementsX, numberElementsY);
xx=nodeCoordinates(:,1);   yy=nodeCoordinates(:,2);
figure
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off

numberNodes=size(xx,1);    % number of nodes

% GDof: global number of degrees of freedom
GDof=5*numberNodes; 

% stiffness and mass matrices
stiffness=formStiffnessMatrixMindlinQ45laminated5dof...
    (GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,CMembranaMembrana,CMembranaFlexao0,CFlexao0Flexao0,CCorte0Corte0);

geometric=formGeometricStiffnessMindlinQ45dof(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness);

% boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=EssentialBC5dof('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);

% buckling analysis

    % perform eigenproblem
    [V1,D1] = eig(stiffness(activeDof,activeDof),geometric(activeDof,activeDof)); 
    D1 = diag(D1);
    % drawing eigenmodes
    numberOfModes=12;
    % sort out eigenvalues
    [D1,ii] = sort(D1); ii = ii(1:numberOfModes); VV = V1(:,ii);
    activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);NNN=size(activeDofW);
    
    % normalize results
    disp('normalized buckling load')
    D1(1:6)'*L*L/e2/h^3
    
    % drawing eigenmodes
    drawingEigenmodes(numberNodes,numberOfModes,NNN,...
        numberElementsX,numberElementsY,...
        L,D1,VV,activeDofW)

end % main function
       
% .............................................................             
    
%   drawing eigenmodes
    
function drawingEigenmodes(numberNodes,numberOfModes,NNN,numx,numy,...
    L,D1,VV,activeDofW)
    
    VVV(1:numberNodes,1:numberOfModes)=0;
    for i=1:numberOfModes
        VVV(activeDofW,i)=VV(1:NNN,i);
    end
%   
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numx+1);
y=linspace(-L,L,numy+1);
% ...............................................
figure
  [xx,yy] = meshgrid(x,y);
  fine = -1:.02:1; 
  [xxx,yyy] = meshgrid(fine,fine);
  uu = zeros(NN,NN);
  [ay,ax] = meshgrid([.56 .04],[.1 .5]); 
  for i = 1:4
    uu = reshape(VVV(1:NN,i),N,N);
    uu = uu/norm(uu(:),inf);
    uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
    subplot('position',[ax(i) ay(i) .38 .38])
    contour(fine,fine,uuu,-0.9:.2:.9)
    colormap(1e-6*[1 1 1]); axis square
    title(['eig = ' num2str(D1(i),'%18.12f')])
  end
%
      figure
  [xx,yy] = meshgrid(x,y);
  fine = -1:.02:1; 
  [xxx,yyy] = meshgrid(fine,fine);
  uu = zeros(NN,NN);
  [ay,ax] = meshgrid([.56 .04],[.1 .5]); 
  for i = 5:8
    uu = reshape(VVV(1:NN,i),N,N);
    uu = uu/norm(uu(:),inf);
    uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
    subplot('position',[ax(i-4) ay(i-4) .38 .38])
    contour(fine,fine,uuu,-0.9:.2:.9)
    colormap(1e-6*[1 1 1]); axis square
    title(['eig = ' num2str(D1(i),'%18.12f')])
  end
%     
  figure
  [xx,yy] = meshgrid(x,y);
  fine = -1:.02:1; 
  [xxx,yyy] = meshgrid(fine,fine);
  uu = zeros(NN,NN);
  [ay,ax] = meshgrid([.56 .04],[.1 .5]); 
  for i = 9:12
    uu = reshape(VVV(1:NN,i),N,N);
    uu = uu/norm(uu(:),inf);
    uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
    subplot('position',[ax(i-8) ay(i-8) .38 .38])
    contour(fine,fine,uuu,-0.9:.2:.9)
    colormap(1e-6*[1 1 1]); axis square
    title(['eig = ' num2str(D1(i),'%18.12f')])
  end
%  
figure  
  % Reshape them to 2D grid, interpolate to finer grid, and plot:
  [xx,yy] = meshgrid(x,y);
  fine_x = -1:.02:1;fine_y = -1:.02:1;  [xxx,yyy] = meshgrid(fine_x,fine_y);
  uu = zeros(NN,NN);
  [ay,ax] = meshgrid([.56 .04],[.1 .5]); clf
  for i = 1:4
    uu = reshape(VVV(1:NN,i),N,N);
    uu = uu/norm(uu(:),inf);
    uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
   subplot(2,2,i)
    meshc(xxx,yyy,uuu)
    title(['eig = ' num2str(D1(i),'%18.12f')])
  end
  end % function drawing eigenmodes
