function [J_mat,N_reduced_x_y_cols]=Jacobian(nodeCoordinates,N_reduced_xi_eta_cols)
% JacobianMatrix    : Jacobian matrix
% XYDerivatives  : derivatives w.r.t. x and y
% naturalDerivatives  : derivatives w.r.t. xi and eta
% nodeCoordinates  : nodal coordinates at element level

J_mat=nodeCoordinates.'*N_reduced_xi_eta_cols;
N_reduced_x_y_cols=N_reduced_xi_eta_cols/J_mat;
