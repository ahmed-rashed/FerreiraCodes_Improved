function [D_vec,F_vec]=solution(prescribedDof,K_assembly,D_vec,F_vec,F_eq_vec)
if nargin<5
    F_eq_vec=zeros(size(F_vec));
end

GDof=length(D_vec);

freeDof=setdiff(1:GDof,prescribedDof);

D_vec(freeDof)=K_assembly(freeDof,freeDof)\(F_vec(freeDof)+F_eq_vec(freeDof));

%nonZeroDof=find(D_vec~=0);
nonZeroDof=union(freeDof,prescribedDof(D_vec(prescribedDof)~=0));

F_vec(prescribedDof)=K_assembly(prescribedDof,nonZeroDof)*D_vec(nonZeroDof)-F_eq_vec(prescribedDof);