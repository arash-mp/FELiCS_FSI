function [v,evec] =  MatlabEigs(n_ev,shift,ncv,maxiter,tol,adjoint_bool)

opts.issym=0;
opts.isreal=0;
opts.tol=tol;

load('A.mat');
load('B.mat');
if adjoint_bool
	[evec,evalue]=eigs(A',B',n_ev,shift,opts);
else
	[evec,evalue]=eigs(A,B,n_ev,shift,opts);
end

for i=1:size(evalue)
    v(i)=evalue(i,i);
end
end
