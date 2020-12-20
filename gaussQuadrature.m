function [w_col,locations_cols]=gaussQuadrature(N)
% Gauss 2D quadrature
% N: Number of integration points

switch N
    case 2
        locations_cols=1/sqrt(3)*[ -1 -1
                                    1 -1
                                    1  1
                                   -1  1];
        w_col=[ 1
                1
                1
                1]; 
    case 1
        locations_cols=[0 0];
        w_col=4;
end