function [m,n,ts]=sum_density_over_bed(d,bd)
%function [m,n,ts]=sum_density_over_bed(d,bd)
%
%IN: d is a Map from chromosomes to alignment density vectors of 
%      50bp bins of alignments
%    bd is a Map from chromosomes to a kX2 matrix of feature start/stops
%
%OUT:m is the total alignment density over the bed intervals
%    n is the number of windows overlapping bed intervals
%    ts is the total number of tags in d

chrs=intersect(d.keys,bd.keys);
m=0;n=0;ts=0;
for i=1:length(chrs)
    ts=ts+sum(d(chrs{i}));
    [mt,nt,~]=sum_on_intervs(d(chrs{i}),bd(chrs{i}),1000);
    m=m+mt;
    n=n+nt;
end