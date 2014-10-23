function w=fuse(rS)
%  function w=fuse(rS)
%
%In: rS is a nxm matrix of m IP density col vecs
%Out: w is the set of weight vectors, ie. the consensus IP is rS*w

if ~matlabpool('size'), matlabpool; end
R=cov(rS,1);
[w,~]=eigs(R,[],1,'lm');
w=w/sum(w);
if any(w<0),w=-w;end
if any(w<0)
    r=rand(size(rS,2),1);
    opts=optimset('Display','off','Algorithm','active-set','UseParallel','always');
    w=fmincon(@(x)-(x'*R*x),r,[],[],[],[],zeros(size(rS,2),1),[],@(x) cons_fun(x),opts);
end
w(w<0)=0;
if any(isnan(w))|any(isinf(w)),w=ones(size(w))/length(w);
else, w=w/sum(w);end
matlabpool close