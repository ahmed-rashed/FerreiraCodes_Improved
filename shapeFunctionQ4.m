function [N_col,N_diff_xi_eta_cols]=shapeFunctionQ4(xi,eta)
% shape function and derivatives for Q4 elements
%
% N_col : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta 
% xi, eta: natural coordinates (-1 ... +1)

N_col=[ (1-xi)*(1-eta)
        (1+xi)*(1-eta)
        (1+xi)*(1+eta)
        (1-xi)*(1+eta)]/4;

N_diff_xi_eta_cols=[-(1-eta)    -(1-xi)
                        1-eta      -(1+xi)
                        1+eta       1+xi
                       -(1+eta)     1-xi]/4;