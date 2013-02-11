function [h p ci zval]=bin_ent_stat(p,q)
%function [h p ci zval]=bin_ent_stat(p,q)
%
%In: p and q are probabilities for Bernoulli trials
%
%Out: h tests the hypothesis that the JS divergence 
%  between p and q is zero, h==1 means the null hypo
%  can be rejected at confidence alpha
%  p is the pvalue, ci is a confidence interval, zval
%  is the z score, help ztest for details

S=cov([p,q;(1-p),(1-q)]);
m=0.5*(p+q);
d=sqrt(0.5*bin_rel_entropy(p,m)+0.5*bin_rel_entropy(q,m));
df=[0;0];
df(1)=(1/(4*sqrt(bin_rel_entropy(p,q))))*log((2*p)/(p+q));
df(2)=(1/(4*sqrt(bin_rel_entropy(p,q))))*log((2*q)/(p+q));
[h,p,ci,zval]=ztest(d,0,sqrt(df'*S*df),.05,'right');
