function H=bin_rel_entropy(p,q)
%  function H=bin_rel_entropy(p,q)
%
%In: p,q are scalors s.t. 0<=p,q<=1
%
%Out: binary relative entropy between p and q

H1=p*log2(p/q);H2=(1-p)*log2((1-p)/(1-q));
if (p==0),H1=0;end
if (p==1),H2=0;end 
H=max(0,H1+H2);