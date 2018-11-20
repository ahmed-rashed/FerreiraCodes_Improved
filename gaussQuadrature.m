% ............................................................. 

    function [weights,locations]=gaussQuadrature(option)
    % Gauss quadrature for Q4 elements
    % option 'complete' (2x2)
    % option 'reduced'  (1x1)
    % locations: Gauss point locations
    % weights: Gauss point weights
        
    switch option
        case 'complete'
    
        locations=...
          [ -0.577350269189626 -0.577350269189626;
             0.577350269189626 -0.577350269189626;
             0.577350269189626  0.577350269189626;
            -0.577350269189626  0.577350269189626];
        weights=[ 1;1;1;1]; 
    
        case 'reduced'
        
        locations=[0 0];
        weights=[4];
    end

    end  % end function gaussQuadrature
    