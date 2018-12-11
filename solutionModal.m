%Created by Ahmed Rashed
function [D_modeShape_cols,w_n_vec]=solutionModal(prescribed_Dofs,D_prescribed,K_assembly,M_assembly,N_modes)
GDof=size(K_assembly,1);
freeDof=setdiff(1:GDof,prescribed_Dofs);

D_modeShape_cols=nan(GDof,N_modes);
D_modeShape_cols(prescribed_Dofs,:)=D_prescribed*ones(1,N_modes);

[D_modeShape_cols(freeDof,:),lambda_mat]=eigs(K_assembly(freeDof,freeDof),M_assembly(freeDof,freeDof),N_modes,'smallestabs','IsSymmetricDefinite',true);
w_n_vec=sqrt(diag(lambda_mat));