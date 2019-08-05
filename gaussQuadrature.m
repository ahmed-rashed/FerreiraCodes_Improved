function [w_col,locations_cols]=gaussQuadrature(option)
% Gauss quadrature for Q4 elements
% option 'complete' (2x2)
% option 'reduced'  (1x1)
% locations: Gauss point locations
% weights: Gauss point weights

switch option
    case 'complete'
        locations_cols=1/sqrt(3)*[ -1 -1
                                    1 -1
                                    1  1
                                   -1  1];
        w_col=[ 1
                1
                1
                1]; 
    case 'reduced'
        locations_cols=[0 0];
        w_col=4;
end