function [J_mat,N_diff_x_y_rows]=Jacobian(nodeCoordinates_rows,N_diff_xi_eta_rows)
% N_diff_xi_eta_rows  : N derivatives w.r.t. xi and eta
% nodeCoordinates  : nodal coordinates at element level

J_mat=N_diff_xi_eta_rows*nodeCoordinates_rows;
N_diff_x_y_rows=J_mat\N_diff_xi_eta_rows;
