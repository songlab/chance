function comp_chance_fdr(data_dir,input_str)
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
        r=strfind(M(tidx(j),:),'Rep'); %which have Rep in ID
        if ~isempty(t)&&~isempty(r)
            S{k}=M(tidx(j),1:r-1); 
            k=k+1;
            cidx=[cidx;j];
        end
    end
    cidx=tidx(cidx);
    [T,~,uidx]=unique(S);
    cidx2={};k=1;
    for j=1:length(T)
        t=find(uidx==j);
        if length(t)>1
            cidx2{k}=cidx(t); 
            k=k+1;
        end
    end
    cidx=cidx2; %cidx : index of control samples with replicates
    if isempty(cidx), continue; end
    for j=1:length(cidx), tidx=setdiff(tidx,cidx{j}); end %tidx : index of treatment samples
    INPUT={}; %INPUT{i} is a m-by-n matrix of n Input replicate samples
    for j=1:length(cidx)
        INPUT{j}=[];
        r=cidx{j};
        for s=1:length(r)
            load(M(r(s),:));
            kz=sample_data.keys;
            d=sample_data(kz{1});
            kz=d.dens.keys; t=[];
            for k=1:length(kz), t=[t;d.dens(kz{k})]; end
            INPUT{j}=[INPUT{j},t];
        end
     end
     IP=[]; %IP is a m-by-k matrix of k IP samples
     for j=1:length(tidx)
         load(M(tidx(j),:)
         kz=sample_data.keys;
         d=sample_data(kz{1});
         kz=d.dens.keys; t=[];
         for k=1:length(kz), t=[t;d.dens(kz{k})]; end
         IP=[IP,t];
     end
end

function [cinull, zvalnull, cialt, zvalalt, pnull, palt]=comp_pvals(IP,INPUT) 
if ~matlabpool('size'), matlabpool; end
nip=size(IP,2); nin=size(INPUT,2);
pnull=[];palt=[];cinull=[];cialt=[];zvalnull=[];zvalalt=[];
cinull_loc=[];cialt_loc=[];zvalnull_loc=[];zvalalt_loc=[];
pnull_loc=[];palt_loc=[];
rn=1;ra=1;
INPUTS={}; 
for j=1:nin %resample each Input replicate, once for each IP sample
    INPUTS{j}=INPUT(:,j);
    parfor k=2:nip, INPUTS{j}=[INPUTS{j},poissrnd(INPUT(:,j))]; end
end
for r=1:nin-1, for s=r+1:nin %all pairwise comparisons of Input replicates
    INPUT1=INPUTS{r}; INPUT2=INPUTS{s};    
    for j=1:nip, parfor k=1:nip %all pairwise comparisons of replicate samplings
            [p,q,~,pnull_loc(rn),~,~,~]=extract_sig(INPUT1(:,j),INPUT2(:,k),[],[]);rn=rn+1;
            [~,~,cinull_loc(:,rn),zvalnull_loc(rn)]=bin_ent_stat(p,q);
            [p,q,~,pnull_loc(rn),~,~,~]=extract_sig(INPUT2(:,k),INPUT1(:,j),[],[]);rn=rn+1;
            [~,~,cinull_loc(:,rn),zvalnull_loc(rn)]=bin_ent_stat(p,q);
            cinull=[cinull cinull_loc];zvalnull=[zvalnull zvalnull_loc];
            pnull=[pnull pnull_loc];
    end, end
end, end
for r=1:nip, for s=1:nin %all pairwise IP-to-Input replicate comparisons
    IP1=IP{r}; INPUT1=INPUTS{s};    
    for j=1:nip, parfor k=1:nip %all pairwise comps of IP to Input replicate samplings
            [p,q,~,palt_loc(ra),~,~,err]=extract_sig(IP1(:,j),INPUT1(:,k),[],[]);
            if any(err==3)
                palt_loc(ra)=1e-30;
                [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,p);
            else
                [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,q);
            end
            ra=ra+1;
            cialt=[cialt cialt_loc];zvalalt=[zvalalt zvalalt_loc];
            palt=[palt palt_loc];
    end, end
end, end

