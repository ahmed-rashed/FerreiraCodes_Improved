%................................................................

% MATLAB codes for Finite Element Analysis
% problem19Buckling.m
% performs buckling analysis
% of Mindlin plates using 3 degrees of freedom per node
% antonio ferreira 2008

clearvars
colordef white

% L: side lenght
L  = 1;     

thickness=0.1;
I=thickness^3/12;h=thickness;

% kapa: shear correction factor
kapa=5/6;         

sigmaX=1/thickness;
sigmaXY=0;
sigmaY=0;
sigmaMatrix=[ sigmaX sigmaXY; sigmaXY sigmaY];

% material properties
syms phi

% reddy material
e2=1;e1=25*e2;g23=0.2*e2;g13=0.5*e2;g12=g13;
miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;

% liew material
e2=1;e1=40*e2;g23=0.5*e2;g13=0.6*e2;g12=g13;
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
    elementNodes,numberNodes,nodeCoordinates,CCorte0Corte0,...
    CFlexao0Flexao0,thickness,I);

[geometric]=...
    formGeometricStiffnessMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,sigmaMatrix,thickness);

% Essential boundary conditions    
    [prescribedDof,activeDof,fixedNodeW]=...
        EssentialBC('ssss',GDof,xx,yy,nodeCoordinates,numberNodes);
  
% buckling analysis ...

    % perform eigenproblem
    [V1,D1] = eig(stiffness(activeDof,activeDof),...
        geometric(activeDof,activeDof)); 
    D1 = diag(D1);
    % drawing eigenmodes
    numberOfModes=12;
    % sort out eigenvalues
    [D1,ii] = sort(D1); ii = ii(1:numberOfModes); 
    VV = V1(:,ii);
    activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
    NNN=size(activeDofW);
    
    % normalize results
    disp('normalized buckling load')
    D1(1:6)'*(L*L/e2/h^3)
    
% drawing eigenmodes
    drawingEigenmodes(numberNodes,numberOfModes,NNN,...
        numberElementsX,numberElementsY,...
        L,D1,VV,activeDofW,CFlexao0Flexao0);
