function EF=forcesInElementGrid(N_elements,elementNodes,xx,yy,E,I,G,J,D_col)

% forces in elements
EF=zeros(6,N_elements);

for iElement=1:N_elements
    iNodes=elementNodes(iElement,:);
    elementDof=[(iNodes(1)-1)*3+1 (iNodes(1)-1)*3+2 (iNodes(1)-1)*3+3 (iNodes(2)-1)*3+1 (iNodes(2)-1)*3+2 (iNodes(2)-1)*3+3] ;
    xa=xx(iNodes(2))-xx(iNodes(1));
    ya=yy(iNodes(2))-yy(iNodes(1));  
    L=sqrt(xa*xa+ya*ya);
    C=xa/L;
    S=ya/L;

    a1 = 12*E*I/(L*L*L);    
    a2 = 6*E*I/(L*L);
    a3 = G*J/L; 
    a4 = 4*E*I/L;
    a5 = 2*E*I/L;

    % stiffness in local axes
    k = [a1 0 a2 -a1 0 a2 ; 0 a3 0 0 -a3 0
         a2 0 a4 -a2 0 a5 ; -a1 0 -a2 a1 0 -a2
         0 -a3 0 0 a3 0; a2 0 a5 -a2 0 a4];

    % transformation matrix
    T_node=[1 0 0
            0 C S
            0 -S C];

    T_element=[T_node zeros(3);zeros(3) T_node];

    % forces in element
    EF(:,iElement)= k*T_element* D_col(elementDof);
end 
