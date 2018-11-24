% Timoshenko beam in bending
% antonio ferreira 2008
%Modified by Ahmed Rashed

clearvars

E=10e7;
nu=0.3;
rho=1;

L=1;
b=1;
h=0.001;
I=b*h^3/12;
kapa=5/6;

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

% computation of the system stiffness matrix 
[K_Assembly,F_equiv]=formStiffnessMassTimoshenkoBeam(GDof,elementNodes,N_nodes,x_col,C,P,rho,I,h);

%force vector
F_col=nan(GDof,1);
F_col(2:N_nodes-1)=0;

% boundary conditions
prescribedDof={ [1 N_nodes N_nodes+[1 N_nodes]]     % clamped-clamped
                [1 N_nodes]            % simply supported-simply supported
                [1 N_nodes+1]};               % clamped at x=0
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
end

% % output D_col/reactions
% outputDisplacementsReactions(D_col,K_Assembly,GDof,prescribedDof)

% drawing deformed shape
subplot(2,1,1)
plot(x_col,D_cols(1:N_nodes,:))
legend({'C-C','S-S','C-F'})