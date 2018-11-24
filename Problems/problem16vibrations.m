% Timoshenko beam in free vibrations
% antonio ferreira 2008
%Modified by Ahmed Rashed

clc
close all
clearvars

E=10e7;
nu=0.30;
rho=1;

L=1;
b=1;
h=0.001;
I=b*h^3/12;
kapa=5/6;
A=b*h;

% constitutive matrix
G=E/2/(1+nu);
C=[E*I   0;
    0    kapa*h*G];

% mesh
N_nodes=101;
x_col=linspace(0,L,N_nodes).';

N_elements=N_nodes-1;
elementNodes=nan(N_elements,2);
elementNodes(:,1)=1:N_elements;
elementNodes(:,2)=2:N_elements+1;

% distributed load
P=-1;

% GDof: global number of degrees of freedom
GDof=2*N_nodes; 

% computation of the system stiffness, force, mass
[K_Assembly,F_equiv,M_Assembly]=formStiffnessMassTimoshenkoBeam(GDof,elementNodes,N_nodes,x_col,C,P,rho,I,h);

%force vector
F_col=nan(GDof,1);
F_col(2:N_nodes-1)=0;

% boundary conditions
prescribedDof={ [1 N_nodes N_nodes+[1 N_nodes]]     % clamped-clamped
                [1 N_nodes]            % simply supported-simply supported
                [1 N_nodes+1]};               % clamped at x=0
titleText={'Clambed-Clambed','Simply-supported-Simply-supported','Clambed-Free'};
            
N_problems=length(prescribedDof);
D_cols=nan(GDof,N_problems);
F_cols=nan(GDof,N_problems);
for ii=1:N_problems
    %displacement vector
    D_cols(prescribedDof{ii},ii)=0;

    %force vector
    freeDOF=setdiff(1:GDof,prescribedDof{ii});
    F_cols(freeDOF,ii)=0;

    % solution
    [D_cols(:,ii),F_cols(:,ii)]=solution(prescribedDof{ii},K_Assembly,D_cols(:,ii),F_cols(:,ii),F_equiv);
    
    % max displacement
    disp(' max displacement')
    min(D_cols(1:N_nodes,ii))
    
    %Normal Modes Analysis
    activeDof=setdiff(1:GDof,prescribedDof{ii});

    [modeShapes,lambda]=eig(K_Assembly(activeDof,activeDof),M_Assembly(activeDof,activeDof));
    w_n=sqrt(lambda);

    N_modes=4;
    v_ModeShape=nan(N_nodes,1);
    [v_activeDof,ind]=setdiff(activeDof,N_nodes+1:GDof);
    figure
    for nn=1:N_modes
        v_ModeShape(setdiff(prescribedDof{ii},N_nodes+1:GDof))=0;
        v_ModeShape(v_activeDof)=modeShapes(ind,nn);
        subplot(N_modes,1,nn)
        plot(x_col,v_ModeShape)
        grid
        ylabel(['mode ',int2str(nn)])
    end
    subplot(N_modes,1,1)
    title([titleText{ii},' mode shapes'])
end

% % output D_col/reactions
% outputDisplacementsReactions(D_col,K_Assembly,GDof,prescribedDof)
