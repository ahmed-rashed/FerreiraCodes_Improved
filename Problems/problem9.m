% Beam (bending only)
% antonio ferreira 2008
% Modified by Ahmed Rashed

close all
clearvars

E=1;

I=1;

% generation of coordinates and connectivities
L=1;
N_Nodes=81;
x_col=linspace(0,L,N_Nodes).';

N_elements=N_Nodes-1;
elementNodes=[(1:N_elements).',(2:N_elements+1).'];

% distributed load
P=-1;
% P=0;

% GDof: global number of degrees of freedom
GDof=2*N_Nodes;

% stiffess matrix and force vector
[K_Assembly,F_equiv]=formStiffnessBernoulliBeam(GDof,N_elements,elementNodes,x_col,E,I,P);

%force vector
F_col=nan(GDof,1);
F_col(3:end-2)=0;
% F_col(39*2+1)=-1;

% boundary conditions
prescribedDof={ [1 2 GDof-1 GDof]     % clamped-clamped
                [1 GDof-1]            % simply supported-simply supported
                [1 2]};               % clamped at x=0
N_problems=length(prescribedDof);
D_cols=nan(GDof,N_problems);
F_cols=nan(GDof,N_problems);
for ii=1:N_problems
    %displacement vector
    D_cols(prescribedDof{ii},ii)=0;
    F_cols(:,ii)=F_col;
    F_cols(setdiff(prescribedDof{1},prescribedDof{ii}),ii)=0;

    % solution
    [D_cols(:,ii),F_cols(:,ii)]=solution(prescribedDof{ii},K_Assembly,D_cols(:,ii),F_cols(:,ii),F_equiv);
end

% % output D_col/reactions
% outputDisplacementsReactions(D_col,stiffness,GDof,prescribedDof)

% drawing deformed shape
subplot(2,1,1)
plot(x_col,D_cols(1:2:end,:))
legend({'C-C','S-S','C-F'})