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
    cidx=cidx2; %index of control samples with replicates
    if isempty(cidx), continue; end
    for j=1:length(cidx), tidx=setdiff(tidx,cidx{j}); end %index of treatment samples
                                                    
end
