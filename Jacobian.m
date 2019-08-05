function [J_mat,N_diff_x_y_cols]=Jacobian(nodeCoordinates,N_diff_xi_eta_cols)
% N_diff_xi_eta_cols  : N derivatives w.r.t. xi and eta
% nodeCoordinates  : nodal coordinates at element level

J_mat=nodeCoordinates.'*N_diff_xi_eta_cols;
N_diff_x_y_cols=N_diff_xi_eta_cols/J_mat;
