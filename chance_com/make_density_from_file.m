function  [d,chr_lens]=make_density_from_file(fname,chr_lens,bin,type)
%function [d,chr_lens]=make_density_from_file(fname,chr_lens,bin,type)
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
%     chr_lens is a map from chromosome ids to chromosome lengths read from
%     the sam/bam file header or from input if not a sam/bam file type

stp=0;k=1;cncl=0;
chunk=1e6; %the number of lines to read in at a time for dsv files
d=[];
try
if strcmp(type, 'bam')|strcmp(type,'sam')
    if ~isdeployed
        javaaddpath('sam-1.64.jar');
        javaaddpath('custombam.jar');
    end
    %get chromosome length information from bam file
    import net.sf.samtools.*
    f=SAMFileReader(java.io.File(fname));
    hd=f.getFileHeader;
    sd=hd.getSequenceDictionary;
    seqs=sd.getSequences;seqs=seqs.toArray;
    if isempty(seqs),disp('This SAM/BAM file is missing a header or the header is corrput.'),return;end
    clear chr_lens
    chr_lens=containers.Map;
    for i=1:length(seqs),chr_lens(char(seqs(i).getSequenceName))=seqs(i).getSequenceLength;end
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
    nreads = double(allData(1));
else 
    f=fopen(fname);
    if strcmp(type,'bed')
        D=textscan(f,'%s%n%*[^\n]',1,'Delimiter','\t');
        d=containers.Map;
        nreads=1;
        while ~feof(f)
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
    elseif strcmp(type,'tagAlign')
        D=textscan(f,'%s%n%*n%s%*[^\n]',1,'Delimiter','\t');
        d=containers.Map;
        seqlen=length(D{3}{:});
        nreads=1;
        while ~feof(f)
            D=textscan(f,['%s%n%*n%' num2str(seqlen) 'c%*[^\n]'],chunk,'Delimiter','\t');
            %D=={chrom,start,seq,phred}
            nreads=nreads+size(D{3},1);
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
    elseif strcmp(type,'bowtie')
        D=textscan(f,'%*s%c%s%n%s%s%*[^\n]',1,'Delimiter','\t');
        d=containers.Map;%create a vector to hold read counts for each chromosome id found
        seqlen=length(D{4}{:});
        nreads=0;
        while ~feof(f)
            D=textscan(f,['%*s%c%s%n%' num2str(seqlen) 'c%' num2str(seqlen) 'c%*[^\n]'],chunk,'Delimiter','\t');
            %D=={strand,chrom,start,seq,phred}
            nreads=nreads+size(D{4},1);
            nidx=(D{1}=='-');%find neg strand reads
            pidx=~nidx;
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
    end
end 
catch me
    disp('error reading file')
    disp(me.message)
    keyboard()
end

