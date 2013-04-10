function  [d,nuc_freq,phred_hist,chr_lens,cncl]=make_density_from_file(fname,chr_lens,bin,type)
%function [d,nuc_freq,phred_hist,chr_lens,cncl]=make_density_from_file(fname,chr_lens,bin,type)
%
%IN:fname is a string holding the file name of the alignments of the form 
%       chr* start *
%   chr_lens is a Map from chromosome id strings to integers giving
%     the chomosome ids and lengths for the given organism
%   bin is the size of each window in bp to bin alignments
%     for now type can be 'bam', s'bed','sam',or 'bowtie', or 'tagAlign'
%
%OUT: d is a map, where d('chr_id')(i) is the number of alignments
%      from file fname which map to the interval
%      ((i-1)*bin+1,i*bin)
%      choose sft==bin for non-overlapping windows
%     nuc_freq is a structure with fields .A .G .C .T, giving vectors of
%     frequencies of nucleotides as a function of sequence index
%     phred_hist is a histogram of phred quality scores
%     chr_lens is a map from chromosome ids to chromosome lengths read from
%     the sam/bam file header or from input if not a sam/bam file type
%     cncl is true if the function returned because the user pressed cancel

stp=0;k=1;cncl=0;
chunk=1e6; %the number of lines to read in at a time
try
if strcmp(type, 'bam')|strcmp(type,'sam')
    if isdeployed
%        javaaddpath(fullfile(ctfroot,'sam-1.64.jar'));
%        javaaddpath(fullfile(ctfroot,'custombam.jar'));
    else
        javaaddpath(fullfile(pwd,'sam-1.64.jar'));
        javaaddpath(fullfile(pwd,'custombam.jar'));
    end
    %get chromosome length information from bam file
    import net.sf.samtools.*
    f=SAMFileReader(java.io.File(fname));
    hd=f.getFileHeader;
    sd=hd.getSequenceDictionary;
    seqs=sd.getSequences;seqs=seqs.toArray;
    clear chr_lens
    chr_lens=containers.Map;
    for i=1:length(seqs),chr_lens(char(seqs(i).getSequenceName))=seqs(i).getSequenceLength;end
    import songlab.*
    allData = CustomBAMMethods.makeDensityFromFile(java.io.File(fname), int32(bin), chr_lens.keys, int32(cell2mat(chr_lens.values')));
    referenceIDsToNames = allData(2);
    i=1:length(referenceIDsToNames);
    referenceKeys = cell(1,length(referenceIDsToNames));
    referenceKeys(i) = referenceIDsToNames(i);
    tempo = allData(3);
    for i=1:length(tempo)
        tempo{i} = double(tempo{i});
    end
    d = containers.Map(referenceKeys,tempo);
    phred_hist = double(allData(5));
    nreads = double(allData(1));
    ATCGvsReadPositionCounts = double(allData(4));
    nuc_freq.A=ATCGvsReadPositionCounts(1,:)/nreads;nuc_freq.T=ATCGvsReadPositionCounts(2,:)/nreads;nuc_freq.G=ATCGvsReadPositionCounts(4,:)/nreads;nuc_freq.C=ATCGvsReadPositionCounts(3,:)/nreads;nuc_freq.N=ATCGvsReadPositionCounts(5,:)/nreads;
    delete(h);
else 
h = waitbar(0,'Progress bar:','Name','Processing reads...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);
dr=dir(fname);
f=fopen(fname);
tmp=textscan(f,'%s',1,'Delimiter','\n');
fclose(f);
lns=dr.bytes/max(length(tmp{1}{:}),100);%estimate the number of lines for a status bar
f=fopen(fname);
if strcmp(type,'bed')
  D=textscan(f,'%s%n%*[^\n]',1,'Delimiter','\t');
  d=containers.Map;
  nreads=1;
  phred_hist=zeros(127,40);
  while ~feof(f)
    cncl=getappdata(h,'canceling');
    if cncl,break,end
    waitbar((k*chunk)/lns);k=k+1;
    D=textscan(f,'%s%n%*[^\n]',chunk,'Delimiter','\t');
    %D=={chrom,start}
    nreads=nreads+size(D{1},1);
    [chrs,~,J]=unique(D{1});%find all the chromosome ids in the file
    %create a vector to hold read counts for each chromosome id found
    for i=1:length(chrs)
        if ~isempty(strfind(lower(chrs{i}),'chr'))&chr_lens.isKey(lower(chrs{i}))
            tmp_d=histc(D{2}(J==i),1:bin:chr_lens(chrs{i}));
            if size(tmp_d,1)<size(tmp_d,2),tmp_d=tmp_d';end
            if d.isKey(chrs{i}), d(chrs{i})=d(chrs{i})+ tmp_d;
            else, d(chrs{i})=tmp_d;end
        end
    end
  end
  nuc_freq.A=zeros(40,1);nuc_freq.C=zeros(40,1);nuc_freq.G=zeros(40,1);nuc_freq.T=zeros(40,1);nuc_freq.N=zeros(40,1);
elseif strcmp(type,'tagAlign')
  D=textscan(f,'%s%n%*n%s%*[^\n]',1,'Delimiter','\t');
  d=containers.Map;
  seqlen=length(D{3}{:});
  nreads=1;
  phred_hist=zeros(127,seqlen);
  a=zeros(1,seqlen);t=a;g=a;c=a;n=a;
  while ~feof(f)
    cncl=getappdata(h,'canceling');
    if cncl,break,end
    waitbar((k*chunk)/lns);k=k+1;
    D=textscan(f,['%s%n%*n%' num2str(seqlen) 'c%*[^\n]'],chunk,'Delimiter','\t');
    %D=={chrom,start,seq,phred}
    nreads=nreads+size(D{3},1);
    a=a+sum(D{3}=='A');t=t+sum(D{3}=='T');
    g=g+sum(D{3}=='G');c=c+sum(D{3}=='C');
    n=n+sum(D{3}=='N'); 
    [chrs,~,J]=unique(D{1});%find all the chromosome ids in the file
    %create a vector to hold read counts for each chromosome id found
    for i=1:length(chrs)
        if ~isempty(strfind(lower(chrs{i}),'chr'))&chr_lens.isKey(lower(chrs{i}))
            tmp_d=histc(D{2}(J==i),1:bin:chr_lens(chrs{i}));
            if size(tmp_d,1)<size(tmp_d,2),tmp_d=tmp_d';end
            if d.isKey(chrs{i}), d(chrs{i})=d(chrs{i})+ tmp_d;
            else, d(chrs{i})=tmp_d;end
        end
    end
  end
  nuc_freq.A=a/nreads;nuc_freq.T=t/nreads;nuc_freq.G=g/nreads;nuc_freq.C=c/nreads;nuc_freq.N=n/nreads;
  %------------depreciated code
  %elseif strcmp(type,'sam')
  %D{1}=textscan(f,'%c%*[^\n]',1);hlns=0;
  %if D{1}=='@',
  %  while D{1}=='@',D=textscan(f,'%c%*[^\n]',1);hlns=hlns+1;end
  %end
  %D=textscan(f,'%*s%n%s%n%*s%*s%*s%*s%*s%s%s%*[^\n]',1,'Delimiter','\t','Headerlines',hlns);
  %d=containers.Map;%create a vector to hold read counts for each chromosome id found
  %seqlen=length(D{4}{:});
  %nreads=0;
  %phred_hist=zeros(127,seqlen);
  %a=zeros(1,seqlen);t=a;g=a;c=a;n=a;
  %while ~feof(f)
  %  cncl=getappdata(h,'canceling');
  %  if cncl,break,end
  %  waitbar((k*chunk)/lns);k=k+1;
  %  D=textscan(f,['%*s%n%s%n%*s%*s%*s%*s%*s%' num2str(seqlen) 'c%' num2str(seqlen) 'c%*[^\n]'],chunk,'Delimiter','\t');
  %  %D=={bit_flg,chrom,start,seq,phred}
  %  nreads=nreads+size(D{4},1);
    %reverse complement all negative strand reads and reverse phredscores
    %  nidx=(bitget(D{1},5)==1);%1 in the fifth bit means - strand,see sam doc
    %  pidx=~nidx;
    %midx=(D{4}=='A');
    %a=a+sum(midx(pidx,:));t=t+fliplr(sum(midx(nidx,:)));
    %midx=(D{4}=='T');
    %t=t+sum(midx(pidx,:));a=a+fliplr(sum(midx(nidx,:)));
    %midx=(D{4}=='G');
    %g=g+sum(midx(pidx,:));c=c+fliplr(sum(midx(nidx,:)));
    %midx=(D{4}=='C');
    %c=c+sum(midx(pidx,:));g=g+fliplr(sum(midx(nidx,:)));
    %midx=(D{4}=='N');
    %n=n+sum(midx(pidx,:))+fliplr(sum(midx(nidx,:)));
    %construct the phred qual score histogram heatmap
    %phtb=double(D{5}); %ASCII value is translated by 33 units (see SAM spec)
    %phred_hist=phred_hist+histc([phtb(pidx,:);fliplr(phtb(nidx,:))],0:126);
    %[chrs,~,J]=unique(D{2});%find all the chromosome ids in the file
    %for i=1:length(chrs)
    %    if ~isempty(strfind(chrs{i},'chr'))
    %   tmp_d=histc(D{3}(J==i),1:bin:(bin-mod(chr_lens(chrs{i}),bin))+chr_lens(chrs{i})+1);
    %        tmp_d = tmp_d(1:length(tmp_d)-1);
    %        if size(tmp_d,1)<size(tmp_d,2),tmp_d=tmp_d';end
    %        if d.isKey(chrs{i}), d(chrs{i})=d(chrs{i})+ tmp_d;
    %        else d(chrs{i})=tmp_d;end
    %    end
        %end
    %end
  %nuc_freq.A=a/nreads;nuc_freq.T=t/nreads;nuc_freq.G=g/nreads;nuc_freq.C=c/nreads;nuc_freq.N=n/nreads;
elseif strcmp(type,'bowtie')
  phtb=[];
  D=textscan(f,'%*s%c%s%n%s%s%*[^\n]',1,'Delimiter','\t');
  d=containers.Map;%create a vector to hold read counts for each chromosome id found
  seqlen=length(D{4}{:});
  nreads=0;
  phred_hist=zeros(127,seqlen);
  a=zeros(1,seqlen);t=a;g=a;c=a;n=a;
  while ~feof(f)
    cncl=getappdata(h,'canceling');
    if cncl,break,end
    waitbar((k*chunk)/lns);k=k+1;
    D=textscan(f,['%*s%c%s%n%' num2str(seqlen) 'c%' num2str(seqlen) 'c%*[^\n]'],chunk,'Delimiter','\t');
    %D=={strand,chrom,start,seq,phred}
    nreads=nreads+size(D{4},1);
    nidx=(D{1}=='-');%find neg strand reads
    pidx=~nidx;
    midx=(D{4}=='A');
    a=a+sum(midx(pidx,:));t=t+fliplr(sum(midx(nidx,:)));
    midx=(D{4}=='T');
    t=t+sum(midx(pidx,:));a=a+fliplr(sum(midx(nidx,:)));
    midx=(D{4}=='G');
    g=g+sum(midx(pidx,:));c=c+fliplr(sum(midx(nidx,:)));
    midx=(D{4}=='C');
    c=c+sum(midx(pidx,:));g=g+fliplr(sum(midx(nidx,:)));
    midx=(D{4}=='N');
    n=n+sum(midx(pidx,:))+fliplr(sum(midx(nidx,:)));
    %construct the phred qual score histogram heatmap
    phtb=double(D{5});
    phred_hist=phred_hist+histc([phtb(pidx,:);fliplr(phtb(nidx,:))],0:126);
    [chrs,~,J]=unique(D{2});%find all the chromosome ids in the file
    for i=1:length(chrs)
        if ~isempty(strfind(lower(chrs{i}),'chr'))&chr_lens.isKey(lower(chrs{i}))
            tmp_d=histc(D{3}(J==i),1:bin:chr_lens(chrs{i}));
            if size(tmp_d,1)<size(tmp_d,2),tmp_d=tmp_d';end
            if d.isKey(chrs{i}),d(chrs{i})=d(chrs{i})+ tmp_d;
            else, d(chrs{i})=tmp_d;end
        end
    end
  end
  nuc_freq.A=a/nreads;nuc_freq.T=t/nreads;nuc_freq.G=g/nreads;nuc_freq.C=c/nreads;nuc_freq.N=n/nreads;
end
if exist('h','var'),delete(h);end
end %close if-else statement
catch me
    disp(me.message)
    for i=1:length(me.stack)
        disp(me.stack.file(i))
        disp(me.stack.name(i))
        disp(me.stack.line(i))
    end
    delete(h);
end

