function [N_reduced_col,N_reduced_xi_eta_cols]=shapeFunctionQ4(xi,eta)
% shape function and derivatives for Q4 elements
%
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta 
% xi, eta: natural coordinates (-1 ... +1)

N_reduced_col=1/4*[ (1-xi)*(1-eta)
            (1+xi)*(1-eta)
            (1+xi)*(1+eta)
            (1-xi)*(1+eta)];

N_reduced_xi_eta_cols=1/4*[-(1-eta), -(1-xi)
                        1-eta,    -(1+xi)
                        1+eta,      1+xi
                        -(1+eta),   1-xi];