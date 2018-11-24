% Beam (bending only)
% antonio ferreira 2008
% Modified by Ahmed Rashed

close all
clearvars

E=1;

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
[K_Assembly,F_equiv]=formStiffnessBernoulliBeam(GDof,elementNodes,x_col,E,I,P);

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

    %force vector
    freeDOF=setdiff(1:GDof,prescribedDof{ii});
    F_cols(freeDOF,ii)=0;

    % solution
    [D_cols(:,ii),F_cols(:,ii)]=solution(prescribedDof{ii},K_Assembly,D_cols(:,ii),F_cols(:,ii),F_equiv);
end

% % output D_col/reactions
% outputDisplacementsReactions(D_col,stiffness,GDof,prescribedDof)

% drawing deformed shape
subplot(2,1,1)
plot(x_col,D_cols(1:2:end,:))
legend({'C-C','S-S','C-F'})