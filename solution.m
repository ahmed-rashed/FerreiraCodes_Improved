function [displacements,force]=solution(prescribedDof,K,displacements,force)

GDof=length(displacements);

freeDof=setdiff(1:GDof,prescribedDof);

displacements(freeDof)=K(freeDof,freeDof)\force(freeDof);

nonZeroDof=find(displacements~=0);

force(prescribedDof)=K(prescribedDof,nonZeroDof)*displacements(nonZeroDof);