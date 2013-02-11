function [m,n,vals] = sum_on_intervs(d,b,res)
%function [m,n,vals] = sum_on_intervs(d,b,res)
%
%In: d is an alignment density
%   b is a nx2 matrix, b(i,:) is the ith interval
%   res is the resolution at which d was generated
%
%Out: m is the sum
%   n is the number of datapoints used to generate it
%   vals are the values of d whose indices map to the intervals of b

b=floor(b/res);m=0;n=0;vals=[];
for i=1:size(b,1)
    v=[];
    if b(i,2)<=length(d)
        v=d(max(1,b(i,1)):b(i,2))';
        vals=[vals,v];
    end
    m=m+sum(v);
    n=n+b(i,2)-b(i,1);
end
end

