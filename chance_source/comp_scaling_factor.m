function [p,q,ht,pval,k,m,sz_ip,sz_input,f,err]=comp_scaling_factor(sig,bak,plot_on)
%[p,q,ht,pval,k,m,sz_ip,sz_input,f,err]=comp_scaling_factor(sig,bak,plot_on)
%
%IN: sig and bak are Maps from chromosome id strings to alignment
%densities for the IP and Input channels respectively
%    plot_on is a bool vector of length 2, plot_on(1) toggles the Lorenz
%    plot and plot_on(2) toggles the linearization plot
%
%Out: 0<=p,q<=1 are percentages, which give the percentage allocation 
%     of tags to background in s1 and s2 respectively
%     ht is 1 if the null hypothesis of the binary entropy divergence test
%     (that the change in entropy between p and q is 0 to linear order) can
%     can be rejected at the %5 significance level, ht ==0 otherwise
%     pval is the pvalue associated with the above test
%     k is the index into the order stat of s1 separating signal bins from background bins
%     m is the total number of bins in the ip
%     sz_ip is the number of reads in the ip
%     sz_input is the number of reads in the input
%     f is a handle to the figure if plot_on, 0 otherwise

f=0;
s=[];b=[];
chrs=intersect(sig.keys,bak.keys);
for i=1:length(chrs)
  s=[s;sig(chrs{i})];
  b=[b;bak(chrs{i})];
end
a1=[];a2=[];
if plot_on(1)
  f1=figure;
  set(f1,'color','w');
  a1=axes;
end
if plot_on(2)
  f2=figure;
  set(f2,'color','w');
  a2=axes;
end
sz_ip=sum(s);sz_input=sum(b);
[p,q,ht,pval,k,m,err]=extract_sig(s,b,a1,a2);