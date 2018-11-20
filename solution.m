function [D_vec,F_vec]=solution(prescribedDof,K,D_vec,F_vec,F_eq_vec)
if nargin<5
    F_eq_vec=zeros(size(F_vec));
end

GDof=length(D_vec);

freeDof=setdiff(1:GDof,prescribedDof);

D_vec(freeDof)=K(freeDof,freeDof)\(F_vec(freeDof)+F_eq_vec(freeDof));

%nonZeroDof=find(D_vec~=0);
nonZeroDof=union(freeDof,prescribedDof(D_vec(prescribedDof)~=0));

F_vec(prescribedDof)=K(prescribedDof,nonZeroDof)*D_vec(nonZeroDof)-F_eq_vec(prescribedDof);