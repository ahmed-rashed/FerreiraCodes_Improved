% 2D problem: thin plate in tension and shear
% antonio ferreira 2008
% Modified by Ahmed Rashed
% This corrects the strange node numbering of Ferreira

clc
clearvars
close all

% materials
E=10e7;
nu=0.3;

rho=1;
h=1;

%Plane stress model
D=E/(1-nu^2)*[  1 nu 0
                nu 1 0
                0 0 (1-nu)/2];

% load
P = 1e6;

%Mesh generation
L_x=5;
L_y=1;
N_elements_X=20;
N_elements_Y=10;
if mod(N_elements_Y,2)~=0,error('N_elements_Y must be even'),end
N_elements=N_elements_X*N_elements_Y;
[nodeCoordinates_rows, elementNodes]=rectangularMesh(L_x,L_y,N_elements_X,N_elements_Y);
N_nodes=size(nodeCoordinates_rows,1);

% GDof: global number of degrees of freedom
GDof=2*N_nodes; 

[K_Assembly,M_Assembly]=formStiffnessMass2D(GDof,elementNodes,nodeCoordinates_rows,D,rho,h);

% boundary conditions 
iNodeLeftEdge=find(nodeCoordinates_rows(:,1)==0);
iNodeMidLeftEdge=iNodeLeftEdge(nodeCoordinates_rows(iNodeLeftEdge,2)==L_y/2);
iNodeRightEdge=find(nodeCoordinates_rows(:,1)==L_x);
iNodeRightTop=iNodeRightEdge(nodeCoordinates_rows(iNodeRightEdge,2)==L_y);
iNodeRightbottom=iNodeRightEdge(nodeCoordinates_rows(iNodeRightEdge,2)==0);

prescribedDof={ [2*iNodeLeftEdge-1;2*iNodeMidLeftEdge]
                [2*iNodeLeftEdge-1;2*iNodeLeftEdge]};

% force vector
F_cols=zeros(GDof,2);

%  Elemental loads (distributed load applied at x=L_x)
F_equiv_cols=zeros(GDof,2);
%problem17.m
F_equiv_cols(2*iNodeRightEdge-1,1)=P*(L_y/N_elements_Y);
F_equiv_cols(2*iNodeRightTop-1,1)=P*(L_y/N_elements_Y)/2;
F_equiv_cols(2*iNodeRightbottom-1,1)=P*(L_y/N_elements_Y)/2;

%problem18.m
F_equiv_cols(2*iNodeRightEdge,2)=P*(L_y/N_elements_Y);
F_equiv_cols(2*iNodeRightTop,2)=P*(L_y/N_elements_Y)/2;
F_equiv_cols(2*iNodeRightbottom,2)=P*(L_y/N_elements_Y)/2;

N_problems=size(F_equiv_cols,2);

% Displacement vector
D_cols=nan(GDof,2);

% Modal analysis
N_modes=4;
D_modeShape_layers=nan(GDof,N_modes,2);
w_n_rows=nan(N_problems,N_modes);

% Drawing
matrixShape=[N_elements_Y+1,N_elements_X+1];
x_mat=reshape(nodeCoordinates_rows(:,1),matrixShape);
y_mat=reshape(nodeCoordinates_rows(:,2),matrixShape);
scaleFactor_vec=[10,0.1];
modeShapescaleFactor=.5;

f1=figure;
for iProblem=1:N_problems
    D_cols(prescribedDof{iProblem},iProblem)=0;
    F_cols(prescribedDof{iProblem},iProblem)=nan;
    
    % Linear static analysis
    [D_cols(:,iProblem),F_cols(:,iProblem)]=solution(prescribedDof{iProblem},K_Assembly,D_cols(:,iProblem),F_cols(:,iProblem),F_equiv_cols(:,iProblem));

    %Normal Modes Analysis
    [D_modeShape_layers(:,:,iProblem),w_n_rows(iProblem,:)]=solutionModal(prescribedDof{iProblem},D_cols(prescribedDof{iProblem},iProblem),K_Assembly,M_Assembly,N_modes);

    %Draw deformed shape
    D_x_mat=reshape(D_cols(1:2:end,iProblem),matrixShape);
    D_y_mat=reshape(D_cols(2:2:end,iProblem),matrixShape);

    x_deformed_mat=x_mat+scaleFactor_vec(iProblem)*D_x_mat;
    y_deformed_mat=y_mat+scaleFactor_vec(iProblem)*D_y_mat;
    D_magnitude_mat=sqrt(D_x_mat.^2+D_y_mat.^2);

    figure(f1);
    subplot(N_problems,1,iProblem)
    h=pcolor(x_deformed_mat,y_deformed_mat,D_x_mat);
    set(h,'EdgeColor',.5*[1,1,1])
    hold on

    mesh(x_mat,y_mat,0*x_mat,'FaceColor','none','EdgeColor','k')
    view(2)
    axis equal
%     axis off
    title(['Deformed shape (scaled ',num2str(scaleFactor_vec(iProblem)),' times)'])
    colorbar
    
    % draw mode shape
    if iProblem==2  %Modes of the 1st problem are difficult to understand
        figure
        for iMode=1:N_modes
            D_x_mat=reshape(D_modeShape_layers(1:2:end,iMode,iProblem),matrixShape);
            D_y_mat=reshape(D_modeShape_layers(2:2:end,iMode,iProblem),matrixShape);
            x_deformed_mat=x_mat+modeShapescaleFactor*D_x_mat;
            y_deformed_mat=y_mat+modeShapescaleFactor*D_y_mat;

            subplot(N_modes,1,iMode)
            mesh(x_mat,y_mat,0*x_mat,'FaceColor','none','EdgeColor',1*[1,0,0])
            hold on
            mesh(x_deformed_mat,y_deformed_mat,0*x_mat,'FaceColor','none','EdgeColor','k')
            view(2)
            axis equal
            axis tight
            ylabel(['$\omega_{',int2str(iMode),'}=',num2str(w_n_rows(iProblem,iMode)),'$'],'interpreter','latex')
        end
        subplot(N_modes,1,1)
        title('Mode shapes')
    end

    % stresses at nodes
    stresses{iProblem}=stresses2D(elementNodes,N_nodes,nodeCoordinates_rows,D_cols(:,iProblem),D_x_mat(:),D_y_mat(:),D,scaleFactor_vec(iProblem));
end
