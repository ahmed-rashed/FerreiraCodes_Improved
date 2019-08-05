% Beam (bending only)
% antonio ferreira 2008
% Modified by Ahmed Rashed

clearvars

E=1;
rho=1;

A=1;
I=1;

% generation of coordinates and connectivities
L=1;
N_nodes=81;
x_col=linspace(0,L,N_nodes).';

N_elements=N_nodes-1;
elementNodes=[(1:N_elements).',(2:N_elements+1).'];

% distributed load
P=-1;
% P=0;

% GDof: global number of degrees of freedom
GDof=2*N_nodes;

% stiffess matrix and force vector
[K_Assembly,F_equiv,M_Assembly]=formStiffnessMassBernoulliBeam(GDof,elementNodes,x_col,E,I,P,rho,A);

% boundary conditions
prescribedDof={ [1 2 GDof-1 GDof]     % clamped-clamped
                [1 GDof-1]            % simply supported-simply supported
                [1 2]};               % clamped at x=0
titleText={'Clamped-Clamped','Simply-supported-Simply-supported','Clambed-Free'};

N_problems=length(prescribedDof);
D_cols=nan(GDof,N_problems);
F_cols=nan(GDof,N_problems);

N_modes=4;
D_modeShape_layers=nan(GDof,N_modes,N_problems);
w_n_rows=nan(N_problems,N_modes);
for iProblem=1:N_problems
    %displacement vector
    D_cols(prescribedDof{iProblem},iProblem)=0;

    %force vector
    freeDOF=setdiff(1:GDof,prescribedDof{iProblem});
    F_cols(freeDOF,iProblem)=0;

    % Linear static analysis
    [D_cols(:,iProblem),F_cols(:,iProblem)]=solution(prescribedDof{iProblem},K_Assembly,D_cols(:,iProblem),F_cols(:,iProblem),F_equiv);
    
    %Normal Modes Analysis
    [D_modeShape_layers(:,:,iProblem),w_n_rows(iProblem,:)]=solutionModal(prescribedDof{iProblem},D_cols(prescribedDof{iProblem},iProblem),K_Assembly,M_Assembly,N_modes);
end

% % output D_col/reactions
% outputDisplacementsReactions(D_col,stiffness,GDof,prescribedDof)

% % draw deformed shape
% figure;
% for iProblem=1:N_problems
%     subplot(N_problems,1,iProblem)
%     plot(x_col,D_cols(1:2:end,iProblem))
%     grid on
%     title([titleText{iProblem},' deformation'])
% end
% 
% % draw mode shape
% for iProblem=1:N_problems
%     figure
%     for iMode=1:N_modes
%         subplot(N_modes,1,iMode)
%         plot(x_col,D_modeShape_layers(1:2:end,iMode,iProblem))
%         grid on
%         ylabel(['$\omega_{',int2str(iMode),'}=',num2str(w_n_rows(iProblem,iMode),'%.4g'),'$'],'interpreter','latex')
%     end
%     subplot(N_modes,1,1)
%     title([titleText{iProblem},' mode shapes'])
% end