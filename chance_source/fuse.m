function w=fuse(rS)
%  function w=fuse(rS)
%
%In: rS is a nxm matrix of m IP density col vecs
%Out: w is the set of weight vectors, ie. the consensus IP is rS*w and the
%consensus Input is rN*w

[w,~]=eigs(rS'*rS,[],1,'lm');
if any(w<0),w=-w;end
if any(w<0)
    r=rand(size(rS,2),1);
    opts=optimset('Display','off','Algorithm','active-set');
    w=fmincon(@(x)-(x'*rS'*rS*x),r,[],[],[],[],zeros(size(rS,2),1),[],[],opts);
end
w=w/sum(w);
