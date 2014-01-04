function [cinull, zvalnull, cialt, zvalalt, pnull, palt]=comp_chance_fdr(data_dir,input_str)
%function comp_chance_fdr(data_dir, input_str)
%
%IN: data_dir - name of directory holding sample BAM files
%    input_str - string identifying control samples
%
%OUT:
%

%dir to find bam files 
%data_dir=/songlab/aaron/research/fun_genom/data/broad_tfbs/normal/
%name of mat file to write containing CHANCE confusion matrix and
%associatd statistics
%fname=broad_tfbs_normal;
%input_str='Control'; %Input-channel identifier, substring of sample file name

cd(data_dir)
pnull=[];palt=[];cinull=[];cialt=[];zvalnull=[];zvalalt=[];
pnull_loc=[];palt_loc=[];cinull_loc=[];cialt_loc=[];zvalnull_loc=[];zvalalt_loc=[];
d=dir('*.mat');
M=strvcat(d.name); %cast cell array as charcter matrix
%identify cell id strings from file name prefix
k=1;
for i=1:size(M,1)
    t=strfind(M(i,:),input_str); %identify control samples
    r=strfind(M(i,:),'Rep'); %which have replicates
    if ~isempty(t)&&~isempty(r)
            cidx(k)=i; 
            cnames{k}=M(i,1:t(1)-1);
            k=k+1;
    end
end
cnames=unique(cnames);
for i=1:length(cnames)
    tidx=strmatch(cnames{i},M); %index samples by cell type
    %index the control samples of cell type cnames{i} which have replicates
    cidx=[];S={};k=1;
    for j=1:length(tidx)
        t=strfind(M(tidx(j),:),input_str); %identify control samples
        r=strfind(M(tidx(j),:),'Rep'); %which have substring 'Rep' in ID
        if ~isempty(t)&&~isempty(r)
            S{k}=M(tidx(j),1:r-1); 
            k=k+1;
            cidx=[cidx;j];
        end
    end
    cidx=tidx(cidx);
    [T,~,uidx]=unique(S);
    cidx2={};k=1;
    for j=1:length(T) %select only those that have reps with same protocol
        t=find(uidx==j);
        if length(t)>1
            cidx2{k}=cidx(t); 
            k=k+1;
        end
    end
    cidx=cidx2; %cidx : index of control samples with replicates
    if isempty(cidx), continue; end
    for j=1:length(cidx), tidx=setdiff(tidx,cidx{j}); end %tidx : index of treatment samples
    INPUT={}; %INPUT{i} is a m-by-n matrix of n Input replicate samples,
              %each index i corresponds to a different protocol
    for j=1:length(cidx)
        INPUT{j}=[];
        r=cidx{j};
        for s=1:length(r) %load all the replicate samples for the ith protocol
            load(deblank(M(r(s),:))); 
            kz=sample_data.keys;
            d=sample_data(kz{1});
            kz=d.dens.keys; t=[];
            for k=1:length(kz), t=[t;d.dens(kz{k})]; end
            try
                INPUT{j}=[INPUT{j},t];
            catch me
                if size(INPUT{j},1)>length(t)
                    t2=zeros(size(INPUT{j},1),1);
                    t2(1:length(t))=t;
                    INPUT{j}=[INPUT{j},t2];
                end
                if size(INPUT{j},1)<length(t)
                    INPUT{j}=[INPUT{j},t(1:size(INPUT{j},1))];
                end
            end
        end
    end
    IP=[]; %IP is a m-by-k matrix of k IP samples
    for j=1:length(tidx) %load all of the IP samples for this cell type
        load(deblank(M(tidx(j),:)));
        if isempty(sample_data), continue; end
        kz=sample_data.keys;
        d=sample_data(kz{1});
        kz=d.dens.keys; t=[];
        for k=1:length(kz), t=[t;d.dens(kz{k})]; end
        try 
            IP=[IP,t];
        catch me
            if size(IP,1)>length(t)
                t2=zeros(size(IP,1),1);
                t2(1:length(t))=t;
                IP=[IP,t2];
            end
            if size(IP,1)<length(t)
                IP=[IP,t(1:size(IP,1))];
            end
        end
    end
    mn=size(IP,1);
    for j=1:length(INPUT)
        t=size(INPUT{j},1);
        if t<mn, mn=t; end
    end
    IP=IP(1:mn,:);
    for j=1:length(INPUT), M=INPUT{j}; INPUT{j}=M(1:mn,:); end
    %resample each Input replicate
    %perform pairwise Input rep-to-Input rep comparisons
    %between all samplings, and perform all pairwise IP-to-Input
    %comparisons. Aggregate the resulting p-values and associated
    %z-scores and confidence intervals
    for j=1:length(INPUT)
        [cinull_loc, zvalnull_loc, cialt_loc, zvalalt_loc, pnull_loc, palt_loc]=comp_pvals(IP,INPUT{j});
        cinull=[cinull cinull_loc]; zvalnull=[zvalnull zvalnull_loc];
        cialt=[cialt cialt_loc]; zvalalt=[zvalalt zvalalt_loc];
        pnull=[pnull pnull_loc]; palt=[palt palt_loc];
    end
end
end %function

function [cinull, zvalnull, cialt, zvalalt, pnull, palt]=comp_pvals(IP,INPUT) 
% function [cinull, zvalnull, cialt, zvalalt, pnull,palt]=comp_pvals(IP,INPUT)
%
%IN: IP - a m-by-n matrix of genome wide IP read densities
%    INPUT - a m-by-k matrix of Input replicate read densities
%
%OUT:  *null - corresponds to Input-to-Input comparisons
%      *alt - corresponds to IP-to-IP comparisons
%      zval* - z-score vector, one entry for each comparison
%      ci* - confidence interval matrix for z-scores
%      p* - p-values for z-scores

if ~matlabpool('size'), matlabpool; end
rp=randperm(size(IP,2));
IP=IP(:,rp(1:min(25,length(rp))));
nip=size(IP,2); nin=size(INPUT,2);
pnull=[];palt=[];cinull=[];cialt=[];zvalnull=[];zvalalt=[];
INPUTS={}; %INPUTS{i} is a matrix of counts, the first column is
           %INPUT(:,i) and the other are resampled from INPUT(:,i)
%resample each Input replicate, so we have the same number of total samples as the IP-to-Input comparison 
nrs=floor((nin*nip-nin*(nin-1)/2)/nin);
for j=1:nin
    M=zeros(size(IP,1),nrs);
    M(:,1)=INPUT(:,j);
    parfor k=2:nrs, M(:,k)=poissrnd(INPUT(:,j)); end
    INPUTS{j}=M;
end
for r=1:nin-1, for s=r+1:nin %all pairwise comparisons of Input replicates
    INPUT1=INPUTS{r}; INPUT2=INPUTS{s};
    n1=size(INPUT1,2); n2=size(INPUT2,2);
    for j=1:n1, parfor k=1:n2 %all pairwise comparisons of replicate samplings
            try
                [p,q,~,pnull_loc,~,~,~]=extract_sig(INPUT1(:,j),INPUT2(:,k),[],[]);
                [~,~,cinull_loc,zvalnull_loc]=bin_ent_stat(p,q);
                [p,q,~,pnull_loc,~,~,~]=extract_sig(INPUT2(:,k),INPUT1(:,j),[],[]);
                [~,~,cinull_loc,zvalnull_loc]=bin_ent_stat(p,q);
                cinull=[cinull cinull_loc]; zvalnull=[zvalnull zvalnull_loc];
                pnull=[pnull pnull_loc];
            catch me, continue; end
    end, end
end, end
for i=1:nin, parfor j=1:nip %all pairwise IP-to-Input comparisons
        try
            [p,q,~,palt_loc,~,~,err]=extract_sig(IP(:,j),INPUT(:,i),[],[]);
            if any(err==3)
                palt_loc=1e-30;
                [~,~,cialt_loc,zvalalt_loc]=bin_ent_stat(p,p);
            else
                [~,~,cialt_loc,zvalalt_loc]=bin_ent_stat(p,q);
            end
            cialt=[cialt cialt_loc]; zvalalt=[zvalalt zvalalt_loc];
            palt=[palt palt_loc];
        catch me
            continue;
        end
end, end
end %function
