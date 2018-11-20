function stresses2Dtruss(numberElements,elementNodes,...
    xx,yy,displacements,E)


% stresses at elements
for e=1:numberElements                          
  indice=elementNodes(e,:);
  elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));
  length_element=sqrt(xa*xa+ya*ya);
  C=xa/length_element;
  S=ya/length_element;   
  sigma(e)=E/length_element* ...
      [-C  -S C S]*displacements(elementDof); 
end    
disp('stresses')
sigma'