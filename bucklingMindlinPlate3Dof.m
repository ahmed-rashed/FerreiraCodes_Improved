function bucklingMindlinPlate3Dof

% this function performs buckling analysis
% of Mindlin plates using 3 degrees of freedom per node

% calls other functions for mesh generation, Gauss quadrature
% shape functions and jacobian matrix

% ......................
clearvars
colordef white

% material properties
% modulusOfElasticity  = Young's modulus
% PoissonRatio         = Poisson's ratio

modulusOfElasticity  = 10920;  % Young
PoissonRatio = 0.30;  % coef. Poisson

% L: side lenght
L  = 1;     

thickness=0.001;
I=thickness^3/12;

% kapa: shear correction factor
kapa=5/6;         

% constitutive matrix
% bending part
C_bending=...
    I*modulusOfElasticity/(1-PoissonRatio^2)*...
    [   1                   PoissonRatio          0 ; 
        PoissonRatio        1                     0 ; 
         0                  0                   (1-PoissonRatio)/2 ];
                 
% shear part
C_shear=...
   kapa*thickness*modulusOfElasticity/2/(1+PoissonRatio)*eye(2);

% initial stress matrix
sigmaX=1/thickness;
sigmaXY=0;
sigmaY=0;
sigmaMatrix=[ sigmaX sigmaXY; sigmaXY sigmaY];
                         
% mesh generation ...
% numberElementsX: number of elements in x
% numberElementsY: number of elements in y
numberElementsX=10;
numberElementsY=10;
% number of elements
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = ...
    rectangularMesh(L, L, numberElementsX, numberElementsY);
xx=nodeCoordinates(:,1);   yy=nodeCoordinates(:,2);
figure
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off

numberNodes=size(xx,1);    % number of nodes
GDof=3*numberNodes;        % total number of DOFs

% stiffness and geometric stiffness matrices
[stiffness]=...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I);

[geometric]=...
    formGeometricStiffnessMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness);

% Essential boundary conditions    
    [prescribedDof,activeDof,fixedNodeW]=...
        EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);
  
% buckling analysis ...

    % perform eigenproblem
    [V1,D1] = eig(stiffness(activeDof,activeDof),geometric(activeDof,activeDof)); 
    D1 = diag(D1);
    % drawing eigenmodes
    numberOfModes=12;
    % sort out eigenvalues
    [D1,ii] = sort(D1); ii = ii(1:numberOfModes); VV = V1(:,ii);
    activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);NNN=size(activeDofW);
    
    % normalize results
    disp('D1(1)/pi/pi/C_bending(1,1)')
    D1(1)/pi/pi/C_bending(1,1)
    D1(1)*pi*pi*C_bending(1,1)
    
    % drawing eigenmodes
    drawingEigenmodes(numberNodes,numberOfModes,NNN,...
        numberElementsX,numberElementsY,...
        L,D1,VV,activeDofW,C_bending)

end % main function
       
% .............................................................             
    
%   drawing eigenmodes
    
function drawingEigenmodes(numberNodes,numberOfModes,NNN,numx,numy,...
    L,D1,VV,activeDofW,C_bending)
    
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
    title(['eig = ' num2str(D1(i)/pi/pi/C_bending(1,1),'%18.12f')])
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
    title(['eig = ' num2str(D1(i)/pi/pi/C_bending(1,1),'%18.12f')])
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
    title(['eig = ' num2str(D1(i)/pi/pi/C_bending(1,1),'%18.12f')])
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
    title(['eig = ' num2str(D1(i)/pi/pi/C_bending(1,1),'%18.12f')])
  end
  end % function drawing eigenmodes
