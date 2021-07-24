function [c,ceq] = constraints(x,max_amp_norm,minEigenvalDist,maxEigenvalDist)

ceq = [];
% we do not impose equality constraints

n = size(x,1);
% twice the number of (real) basis functions
halfn=round(n/2);

tp=[[x(1:halfn,1)+1i*x((halfn+1):end,1)];[x(1:halfn,1)-1i*x((halfn+1):end,1)]];
% all pairs of complex conjugate eigenvalues

tot_mat=abs(repmat(tp,1,n)-repmat(tp.',n,1));
% Matrix of distances between eigenvalues (distances without repeats are in in the indices triu(ones(n),1)==1 )


c = [((sum((x(:,2)).^2))-max_amp_norm.^2) ; (minEigenvalDist-min(tot_mat(triu(ones(n),1)==1))) ; (max(tot_mat(triu(ones(n),1)==1))-maxEigenvalDist)];
% Max square norm of the amplitudes will be max_amp_norm.
% minEigenvalDist is the minimal distance between eigenvalues
% maxEigenvalDist is the maximal distance between eigenvalues
