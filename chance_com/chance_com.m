function out=chance_com(subr,varargin)
%function chance_com(subr,varargin)
%
%IN: subr is a subroutine to execute,
%    varargin: variable argument input of parameter value pairs,
%    see README for usage
%
%OUT:

out=0;
cmds={'binData','IPStrength','multiIPNorm','compENCODE','spectrum'};
if ~ismember(subr,cmds)|isempty(varargin),disp_help();return;end
if strcmp(subr,'binData')
    options = containers.Map({'-p','-b','-t','-s','-o','-f'},{[],[],[],[],[],[]});
elseif strcmp(subr,'IPStrength')
    options = containers.Map({'-p','-b','-t','-o','--ipfile','--ipsample','--inputfile','--inputsample'},{[],[],[],[],[],[],[],[]});
elseif strcmp(subr,'multiIPNorm')
    options = containers.Map({'-p'},{[]});
elseif strcmp(subr,'compENCODE')
    options = containers.Map({'-p','-b','-t','-o','-e','--ipfile','--ipsample','--inputfile','--inputsample'},{[],[],[],[],[],[],[],[],[]});
elseif strcmp(subr,'spectrum')
    options = containers.Map({'-p','-b','-t','-s','-o','-f'},{[],[],[],[],[],[]});
end
optionNames = options.keys;
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2,disp_help();return;end
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1});
   if any(strmatch(inpName,optionNames))
      options(inpName) = pair{2};
   else
      disp([inpName ' is not a recognized parameter name'])
      disp_help();
      return;
   end
end
if strcmp(subr,'binData')
    if ~isempty(options('-p'))
        %parameter file must be comma separated values
        %input_file_name,output_file_name,sample_id,build,file_type
        try, f=fopen(options('-p'));D=textscan(f,'%s%s%s%s%s','Delimiter',',');fclose(f);
        catch me, disp('error opening/parsing parameter file, please check file...'),end
        load('hg18lengths.mat');hg18_chr_lens=chr_lens;
        load('hg19lengths.mat');hg19_chr_lens=chr_lens;
        load('mm9lengths.mat');mm9_chr_lens=chr_lens;
        clear chr_lens;
        for i=1:length(D{1})
            fin{i}=D{1}{i};outf{i}=D{2}{i};smp_id{i}=D{3}{i};
            bld{i}=D{4}{i};lenf{i}=[bld{i},'lengths.mat'];typ{i}=D{5}{i};
            if strcmpi(bld{i},'hg18'),chr_lens{i}=hg18_chr_lens;
            elseif strcmpi(bld{i},'hg19'), chr_lens{i}=hg19_chr_lens;
            else, chr_lens{i}=mm9_chr_lens;end
        end
        smpd=par_bin_data(fin,smp_id,bld,chr_lens,typ);
        for i=1:length(D{1})
            sample_data=smpd{i};save(outf{i},'sample_data');
            smp=sample_data(smp_id{i});
            otf=outf{i};lst=strfind(otf,'.')-1;
            if isempty(lst),lst=length(otf);end
            otf=otf(1:lst)
            csvwrite([otf,'_phred.csv'],smp.phred)
            csvwrite([otf,'_Afreq.csv'],smp.nuc_freq.A)
            csvwrite([otf,'_Gfreq.csv'],smp.nuc_freq.G)
            csvwrite([otf,'_Cfreq.csv'],smp.nuc_freq.C)
            csvwrite([otf,'_Tfreq.csv'],smp.nuc_freq.T)
            csvwrite([otf,'_Nfreq.csv'],smp.nuc_freq.N)
        end
    else
        bld=options('-b');
        if isempty(bld)||~ismember(bld,{'hg18','hg19','mm9','tair10'})
            disp('valid build options are hg18, hg19, mm9, or tair10')
            return;
        end
        typ=options('-t');
        if isempty(typ)||~ismember(typ,{'bam','sam','bed','bowtie','tagAlign','mat'})
            disp('valid file type options are bam, sam, bed, bowtie, tagAlign, or mat')
            return;
        elseif strcmp(typ,'mat'), disp('this is redundant...exiting'),return;end
        load([bld 'lengths.mat']);
        [d,nuc_freq,phred_hist,~,~]=make_density_from_file(options('-f'),chr_lens,1000,typ);
        disp('finished binning reads...')
        chrs=d.keys;n=0;
        for i=1:length(chrs),n=n+sum(d(chrs{i}));end
        smp.nreads=n;smp.genome=bld;
        smp.dens=d;smp.nuc_freq=nuc_freq;smp.phred=phred_hist;
        sample_data(options('-s'))=smp;
        if isKey(options,'-o'), outf=options('-o');
        else, outf='new_sample.mat';end
        save(outf,'sample_data');
        %lst=strfind(outf,'.')-1;
        %if isempty(lst),lst=length(outf);end
        %outf=outf(1:lst)
        %csvwrite([outf,'_phred.csv'],smp.phred)
        %csvwrite([outf,'_Afreq.csv'],smp.nuc_freq.A)
        %csvwrite([outf,'_Gfreq.csv'],smp.nuc_freq.G)
        %csvwrite([outf,'_Cfreq.csv'],smp.nuc_freq.C)
        %csvwrite([outf,'_Tfreq.csv'],smp.nuc_freq.T)
        %csvwrite([outf,'_Nfreq.csv'],smp.nuc_freq.N)
    end
elseif strcmp(subr,'IPStrength')
    if ~isempty(options('-p'))
        %parameter file must be comma separated values
        %IP_file_name,Input_file_name,IP_sample_id,Input_sample_id,output_file,build,file_type
        try, f=fopen(options('-p'));D=textscan(f,'%s%s%s%s%s%s%s','Delimiter',',');fclose(f);
        catch me, disp('error opening/parsing parameter file, please check file...'),end
        load('hg18lengths.mat');hg18_chr_lens=chr_lens;
        load('hg19lengths.mat');hg19_chr_lens=chr_lens;
        load('mm9lengths.mat');mm9_chr_lens=chr_lens;
        clear chr_lens;
        for i=1:length(D{1})
            ipf{i}=D{1}{i};inputf{i}=D{2}{i};ip_smp_id{i}=D{3}{i};input_smp_id{i}=D{4}{i};
            outf{i}=D{5}{i};bld{i}=D{6}{i};typ{i}=D{7}{i};
            if strcmpi(bld{i},'hg18'),chr_lens{i}=hg18_chr_lens;
            elseif strcmpi(bld{i},'hg19'), chr_lens{i}=hg19_chr_lens;
            else, chr_lens{i}=mm9_chr_lens;end
        end
        midx=find(strcmp(typ,'mat'));
        nidx=setdiff([1:length(typ)],midx);
        if ~isempty(nidx)
            input_smpd(nidx)=par_bin_data(inputf(nidx),input_smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
            ip_smpd(nidx)=par_bin_data(ipf(nidx),ip_smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
        end
        for i=1:length(midx)
            sample_data=[];
            load(ipf{midx(i)});
            ip_smpd{midx(i)}=sample_data;
            sample_data=[];
            load(inputf{midx(i)});
            input_smpd{midx(i)}=sample_data;
        end
        batch_ip_strength(input_smpd,ip_smpd,inputf,ipf,input_smp_id,ip_smp_id,outf);
    else
        %valid options
        %'-b','-t','-o','--ipfile','--ipsample','--inputfile','--inputsample'
        bld=options('-b');
        if isempty(bld)||~ismember(bld,{'hg18','hg19','mm9','tair10'})
            disp('valid build options are hg18, hg19, mm9, or tair10')
            return;
        end
        typ=options('-t');
        if isempty(typ)||~ismember(typ,{'bam','sam','bed','bowtie','tagAlign','mat'})
            disp('valid file type options are bam, sam, bed, bowtie, tagAlign, or mat')
            return;
        end
        load([bld 'lengths.mat']);
        if strcmp(options('-t'),'mat')
            load(options('--inputfile'));
            input_sample=sample_data(options('--inputsample'));
            load(options('--ipfile'));
            ip_sample=sample_data(options('--ipsample'));
        else
            %load ip file
            [d,nuc_freq,phred_hist,~,~]=make_density_from_file(options('--ipfile'),chr_lens,1000,typ);
            chrs=d.keys;n=0;
            for i=1:length(chrs),n=n+sum(d(chrs{i}));end
            smp.nreads=n;smp.genome=bld;
            smp.dens=d;smp.nuc_freq=nuc_freq;smp.phred=phred_hist;
            ip_sample=smp;
            %load input file
            [d,nuc_freq,phred_hist,~,~]=make_density_from_file(options('--inputfile'),chr_lens,1000,typ);
            chrs=d.keys;n=0;
            for i=1:length(chrs),n=n+sum(d(chrs{i}));end
            smp.nreads=n;smp.genome=bld;
            smp.dens=d;smp.nuc_freq=nuc_freq;smp.phred=phred_hist;
            input_sample=smp;
        end
        [s,fd,ht,k,m,sz_ip,sz_input,p,q]=ip_strength(input_sample,ip_sample);
        try, f=fopen(options('-o'),'w');catch me, disp(['error opening output file ' options('-o')]),end
        fprintf(f,'IP_file,%s\n',options('--ipfile'));fprintf(f,'Input_file,%s\n',options('--inputfile'));
        fprintf(f,'IP_sample_id,%s\n',options('--ipsample'));fprintf(f,'Input_sample_id,%s\n',options('--inputsample'));
        fprintf(f,'pass,%g\n',ht);
        fprintf(f,'fdr,%g\n',fd('all'));fprintf(f,'tfbs_normal_fdr,%g\n',fd('tfbs_normal'));
        fprintf(f,'histone_normal_fdr,%g\n',fd('histone_normal'));fprintf(f,'tfbs_cancer_fdr,%g\n',fd('tfbs_cancer'));
        fprintf(f,'histone_cancer_fdr,%g\n',fd('histone_cancer'));fprintf(f,'percent_genome_enriched,%g\n',(100-100*k/m));
        fprintf(f,'input_scaling_factor,%g\n',(p*sz_ip)/(q*sz_input));
        fprintf(f,'differential_percentage_enrichment,%g\n',100*(q-p));
        for j=1:length(s),fprintf(f,'%s\n',s{j});,end
        if f~=-1,fclose(f);end
    end
elseif strcmp(subr,'multiIPNorm')
    if ~isempty(options('-p'))
        %parameter file must be comma separated values
        %file_name,sample_id,output_file,build,file_type
        try, f=fopen(options('-p'));D=textscan(f,'%s%s%s%s%s','Delimiter',',');fclose(f);
        catch me, disp('error opening/parsing parameter file, please check file...'),end
        load('hg18lengths.mat');hg18_chr_lens=chr_lens;
        load('hg19lengths.mat');hg19_chr_lens=chr_lens;
        load('mm9lengths.mat');mm9_chr_lens=chr_lens;
        clear chr_lens;
        outf=D{3}{1};
        for i=1:length(D{1})
            ipf{i}=D{1}{i};smp_id{i}=D{2}{i};bld{i}=D{4}{i};typ{i}=D{5}{i};
            if strcmpi(bld{i},'hg18'),chr_lens{i}=hg18_chr_lens;
            elseif strcmpi(bld{i},'hg19'), chr_lens{i}=hg19_chr_lens;
            else, chr_lens{i}=mm9_chr_lens;end
        end
        ip_samples_id=smp_id;
        midx=find(strcmp(typ,'mat'));
        nidx=setdiff([1:length(typ)],midx);
        if ~isempty(nidx)
            smpd(nidx)=par_bin_data(ipf(nidx),smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
        end
        for i=1:length(midx)
            sample_data=[];
            load(ipf{midx(i)});
            smpd{midx(i)}=sample_data;
        end
        num_samples=length(smpd);
        for i=1:num_samples
            tmp=smpd{i};
            ip_samples{i}=tmp(smp_id{i});
            ipt=ip_samples{i};ipt=ipt.dens;
            if i==1,kz=ipt.keys;
            else, kz=intersect(kz,ipt.keys);end
        end
        %create the genome wide read density vector for each sample
        for i=1:num_samples
            ipt=ip_samples{i};ipt=ipt.dens;
            ipl=[];for j=1:length(kz),ipl=[ipl;ipt(kz{j})];end
            rS(:,i)=ipl;
        end
        mcnt=mean(sum(rS));ip_depths=sum(rS);
        for i=1:length(ip_depths),rS(:,i)=rS(:,i)/ip_depths(i);end %normalize first by sequencing depth
        w=fuse(rS);%compute weights via signal combiner
        cons_ip=rS*w;
        s1=cons_ip;
        m=length(s1);
        [ss1,idx]=sort(s1);
        cs1=cumsum(ss1);
        gz=max(find(ss1==0))+1;
        ss1_cut=ss1(gz:end);
        cs1_cut=cumsum(ss1_cut);
        t={};
        txt_out={'Differential enrichment:',''};txt_idx=1;
        t{length(t)+1}=txt_out{end-1};
        load('div_fdr.mat');%loads t==pvals and fdrs
        for i=1:num_samples 
            s2=rS(:,i);%reorder the sample with respect to the consensus order stat
            s2r=s2(idx);
            CS2(:,i)=cumsum(s2r);
            %cut the leading zeros from signal is
            s2r_cut=s2r(gz:end);
            cs2_cut=cumsum(s2r_cut);
            %compute the point of maximal difference for the cut dataset
            [~,k(i)]=max(abs(cs1_cut/cs1_cut(end)-cs2_cut/cs2_cut(end)));
            k(i)=k(i)+gz;
            %compute the sig/bak cutoff, adding the zero bins back in and plot
            p(i)=cs1(k(i))/cs1(end);q(i)=CS2(k(i),i)/CS2(end,i);
            [~,pval(i),~,~]=bin_ent_stat(p(i),q(i));
            if pval(i)>=max(pvs),fd(i)=max(fdrs);
            elseif pval(i)<=min(pvs), fd(i)=min(fdrs);
            else,fd(i)=interp1(pvs,fdrs,pval(i));end
            ip_scale(i)=(mcnt/ip_depths(i))*(p(i)/q(i));
            %text output
            txt_out{txt_idx}=[ip_samples_id{i} ':'];txt_idx=txt_idx+1;t{length(t)+1}=[ip_samples_id{i} ':'];
            txt_out{txt_idx}=['Scale ' ip_samples_id{i} ' by ' num2str(ip_scale(i))];txt_idx=txt_idx+1;
            t{length(t)+1}=txt_out{end};
            if fd(i)<=0.05
                txt_out{txt_idx}='Significant differential enrichment from consensus,';txt_idx=txt_idx+1;
                t{length(t)+1}=txt_out{end};
                txt_out{txt_idx}=['q-value (pFDR): ' num2str(fd(i))];txt_idx=txt_idx+1;
                t{length(t)+1}=txt_out{end};
            else
                txt_out{txt_idx}='Differential enrichment from consensus';txt_idx=txt_idx+1;t{length(t)+1}=txt_out{end};
                txt_out{txt_idx}='does not appear significant,';txt_idx=txt_idx+1;t{length(t)+1}=txt_out{end};
                txt_out{txt_idx}=['q-value (pFDR): ' num2str(fd(i))];txt_idx=txt_idx+1;t{length(t)+1}=txt_out{end};
            end
            txt_out{txt_idx}='';txt_idx=txt_idx+1;
        end
        %compute pairwise differential enrichment between samples at the point of maximum
        %divergence from consensus 
        diff_enrich=zeros(num_samples);diff_genom=zeros(num_samples);
        kdiff=ones(num_samples);pdiff=ones(num_samples);qdiff=ones(num_samples);
        for i=1:num_samples
            for j=i+1:num_samples
                [~,kdiff(i,j)]=max(abs(CS2(:,i)/CS2(end,i)-CS2(:,j)/CS2(end,j)));
                pdiff(i,j)=CS2(kdiff(i,j),i)/CS2(end,i);qdiff(i,j)=CS2(kdiff(i,j),j)/CS2(end,j);
                [~,pvaldiff(i,j),~,~]=bin_ent_stat(pdiff(i,j),qdiff(i,j));    
                if pvaldiff(i,j)>=max(pvs),fddiff(i,j)=max(fdrs);
                elseif pvaldiff(i,j)<=min(pvs), fddiff(i,j)=min(fdrs);
                else,fddiff(i,j)=interp1(pvs,fdrs,pvaldiff(i,j));end
                if fddiff(i,j)<=0.05 %if the differential enrichment is not significant than the percent enrichment figure is not meaningful
                    diff_enrich(i,j)=abs(pdiff(i,j)-qdiff(i,j));diff_enrich(j,i)=diff_enrich(i,j);
                    diff_genom(i,j)=1-kdiff(i,j)/m;diff_genom(j,i)=diff_genom(i,j);
                end
            end
        end
        for i=1:length(ip_samples_id)
            tmp_str=ip_samples_id{i};
            stp=strfind(tmp_str,' - # reads')-1;
            if isempty(stp),stp=min(length(tmp_str),10);end
            ssmp_id{i}=tmp_str(1:stp);
        end
        try, f=fopen(outf,'w');catch me, disp(['error opening output file ' outf]),end
        for j=1:length(t),fprintf(f,'%s\n',t{j});,end
        fprintf(f,'\n');
        fprintf(f,'differential enrichment:\n');
        for i=1:length(ssmp_id)-1,fprintf(f,'%s,',ssmp_id{i});end
        fprintf(f,'%s\n',ssmp_id{end});
        for i=1:size(diff_genom,1)
            for j=1:size(diff_genom,2)-1
                fprintf(f,'%g,',diff_genom(i,j));
            end
            fprintf(f,'%g\n',diff_genom(i,size(diff_genom,2)));
        end
        if f~=-1,fclose(f);end
    end
elseif strcmp(subr,'compENCODE')
    if ~isempty(options('-p'))
        %parameter file must be comma separated values
        %IP_file_name,Input_file_name,IP_sample_id,Input_sample_id,exp_id,output_file,build,file_type
        try, f=fopen(options('-p'));D=textscan(f,'%s%s%s%s%s%s%s%s','Delimiter',',');fclose(f);
        catch me, disp('error opening/parsing parameter file, please check file...'),end
        load('hg18lengths.mat');hg18_chr_lens=chr_lens;
        load('hg19lengths.mat');hg19_chr_lens=chr_lens;
        load('mm9lengths.mat');mm9_chr_lens=chr_lens;
        clear chr_lens;
        for i=1:length(D{1})
            ipf{i}=D{1}{i};inputf{i}=D{2}{i};ip_smp_id{i}=D{3}{i};input_smp_id{i}=D{4}{i};
            exp_id{i}=D{5}{i};outf{i}=D{6}{i};bld{i}=D{7}{i};typ{i}=D{8}{i};
            if strcmpi(bld{i},'hg18'),chr_lens{i}=hg18_chr_lens;
            elseif strcmpi(bld{i},'hg19'), chr_lens{i}=hg19_chr_lens;
            else, chr_lens{i}=mm9_chr_lens;end
        end
        midx=find(strcmp(typ,'mat'));
        nidx=setdiff([1:length(typ)],midx);
        if ~isempty(nidx)
            input_smpd(nidx)=par_bin_data(inputf(nidx),input_smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
            ip_smpd(nidx)=par_bin_data(ipf(nidx),ip_smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
        end
        for i=1:length(midx)
            sample_data=[];
            load(ipf{midx(i)});
            ip_smpd{midx(i)}=sample_data;
            sample_data=[];
            load(inputf{midx(i)});
            input_smpd{midx(i)}=sample_data;
        end
        batch_comp_encode(input_smpd,ip_smpd,inputf,ipf,input_smp_id,ip_smp_id,exp_id,bld,outf) 
    else
        %valid options
        %-b,-t,-o,-e,--ipfile,--ipsample,--inputfile,--inputsample
        bld=options('-b');
        if isempty(bld)||~ismember(bld,{'hg18','hg19','mm9','tair10'})
            disp('valid build options are hg18, hg19, mm9, or tair10')
            return;
        end
        typ=options('-t');
        if isempty(typ)||~ismember(typ,{'bam','sam','bed','bowtie','tagAlign','mat'})
            disp('valid file type options are bam, sam, bed, bowtie, tagAlign, or mat')
            return;
        end
        exp_id=options('--ipsample');
        if isempty(exp_id),exp_id='my_sample';end
        load([bld 'lengths.mat']);
        if strcmp(options('-t'),'mat')
            load(options('--inputfile'));
            input_sample=sample_data(options('--inputsample'));
            load(options('--ipfile'));
            ip_sample=sample_data(options('--ipsample'));
        else
            [ip_dens,nuc_freq,phred_hist,~,~]=make_density_from_file(options('--ipfile'),chr_lens,1000,typ);
            [input_dens,nuc_freq,phred_hist,~,~]=make_density_from_file(options('--inputfile'),chr_lens,1000,typ);
        end
        if strcmp(bld,'hg19'),load('hg19_sn_models.mat');
        else, load('mm9_sn_models.mat');end
        [od,p]=find_tf_binding_odds(options('-e'),ip_dens,input_dens,1000,tf_beds,tf_dists);
        t{1}=['Signal to noise ratio (SNR): ' num2str(od)];
        t{2}=['The probability of observing the given SNR or less in the ENCODE database: ' num2str(p)];
        t{3}=['A small probability indicates your data differs greatly from ENCODE datasets'];
        try, f=fopen(options('-o'),'w');catch me, disp(['error opening output file ' options('-o')]),end
        fprintf(f,'IP_file,%s\n',options('--ipfile'));fprintf(f,'Input_file,%s\n',options('--inputfile'));
        fprintf(f,'IP_sample_id,%s\n',options('--ipsample'));fprintf(f,'Input_sample_id,%s\n',options('--inputsample'));
        fprintf(f,'experiment_id,%s\n',exp_id);fprintf(f,'build,%s\n',bld);
        fprintf(f,'odds_ratio,%g\n',od);fprintf(f,'probability,%g\n',p);
        for j=1:length(t),fprintf(f,'%s\n',t{j});,end
        if f~=-1,fclose(f);end
    end
elseif strcmp(subr,'spectrum')
    if ~isempty(options('-p'))
        %parameter file must be comma separated values
        %file_name,sample_id,output_file,build,file_type
        try, f=fopen(options('-p'));D=textscan(f,'%s%s%s%s%s','Delimiter',',');fclose(f);
        catch me, disp('error opening/parsing parameter file, please check file...'),end
        load('hg18lengths.mat');hg18_chr_lens=chr_lens;
        load('hg19lengths.mat');hg19_chr_lens=chr_lens;
        load('mm9lengths.mat');mm9_chr_lens=chr_lens;
        clear chr_lens;
        for i=1:length(D{1})
            ipf{i}=D{1}{i};smp_id{i}=D{2}{i};
            outf{i}=D{3}{i};bld{i}=D{4}{i};typ{i}=D{5}{i};
            if strcmpi(bld{i},'hg18'),chr_lens{i}=hg18_chr_lens;
            elseif strcmpi(bld{i},'hg19'), chr_lens{i}=hg19_chr_lens;
            else, chr_lens{i}=mm9_chr_lens;end
        end
        midx=find(strcmp(typ,'mat'));
        nidx=setdiff([1:length(typ)],midx);
        if ~isempty(nidx)
            smpd(nidx)=par_bin_data(ipf(nidx),smp_id(nidx),bld(nidx),chr_lens(nidx),typ(nidx));
        end
        for i=1:length(midx)
            sample_data=[];
            load(ipf{midx(i)});
            smpd{midx(i)}=sample_data;
        end
        batch_spectrum(smpd,ipf,smp_id,bld,outf) 
    else
        %valid options
        %-b,-t,-o,-f,-s
        bld=options('-b');
        if isempty(bld)||~ismember(bld,{'hg18','hg19','mm9','tair10'})
            disp('valid build options are hg18, hg19, mm9, or tair10')
            return;
        end
        typ=options('-t');
        if isempty(typ)||~ismember(typ,{'bam','sam','bed','bowtie','tagAlign','mat'})
            disp('valid file type options are bam, sam, bed, bowtie, tagAlign, or mat')
            return;
        end
        load([bld 'lengths.mat']);
        outf=options('-o');
        if isempty(outf),disp('no output file specified'),return;end
        if strcmp(options('-t'),'mat')
            load(options('-f'));
            smp=sample_data(options('-s'));
            dens=smp.dens;
        else
            [dens,nuc_freq,phred_hist,~,~]=make_density_from_file(options('-f'),chr_lens,1000,typ);
        end
        Smpl=[];chrs=dens.keys;
        for j=1:length(chrs),Smpl=[Smpl;dens(chrs{j})];end
        [c,l]=wavedec(Smpl,15,'haar');
        [eau,ed_inp]=wenergy(c,l);
        smp_hist=ed_inp'/sum(ed_inp);
        t{1}=['apx_coef_energy_user,' num2str(eau)];
        d=fitdist(max(Smpl,1),'gamma');
        sim_data=poissrnd(d.random(length(Smpl),1));
        [c,l]=wavedec(sim_data(find(~isnan(sim_data))),15,'haar');
        [eas,ed_sim]=wenergy(c,l);
        sim_hist=ed_sim'/sum(ed_sim);
        t{2}=['apx_coef_energy_sim' num2str(eas)];
        try, f=fopen(outf,'w');catch me, disp(['error opening output file ' outf]),end
        fprintf(f,'user_file,%s\n',options('-f'));fprintf(f,'sample_id,%s\n',options('-s'));
        fprintf(f,'bld,%s\n',options('-b'));
        for j=1:length(t),fprintf(f,'%s\n',t{j});,end
        if f~=-1,fclose(f);end
        otf=outf;
        lst=strfind(otf,'.')-1;
        if isempty(lst),lst=length(otf);end
        otf=otf(1:lst);
        csvwrite([otf,'_user_hist.csv'],smp_hist);
        csvwrite([otf,'_sim_hist.csv'],sim_hist);
    end      
end

function out=batch_spectrum(smpd,ipf,smp_id,bld,outf)
out=0;
s=cell(length(outf),1);
smp_hist=s;sim_hist=s;
if ~matlabpool('size'),matlabpool;end
parfor i=1:length(s)
    t=cell(2,1);
    smp=containers.Map;
    smp=smpd{i};
    sample_data=smp(smp_id{i});
    Smpl=[];dens=sample_data.dens;chrs=dens.keys;
    for j=1:length(chrs),Smpl=[Smpl;dens(chrs{j})];end
    [c,l]=wavedec(Smpl,15,'haar');
    [eau,ed_inp]=wenergy(c,l);
    smp_hist{i}=ed_inp'/sum(ed_inp);
    t{1}=['apx_coef_energy_user,' num2str(eau)];
    d=fitdist(max(Smpl,1),'gamma');
    sim_data=poissrnd(d.random(length(Smpl),1));
    [c,l]=wavedec(sim_data(find(~isnan(sim_data))),15,'haar');
    [eas,ed_sim]=wenergy(c,l);
    sim_hist{i}=ed_sim'/sum(ed_sim);
    t{2}=['apx_coef_energy_sim' num2str(eas)];
    s{i}=t;
end
if matlabpool('size'),matlabpool close;end
for i=1:length(s)
    try, f=fopen(outf{i},'w');catch me, disp(['error opening output file ' outf{i}]),end
    fprintf(f,'file,%s\n',ipf{i});fprintf(f,'sample_id,%s\n',smp_id{i});
    fprintf(f,'build,%s\n',bld{i});
    t=s{i};for j=1:length(t),fprintf(f,'%s\n',t{j});end
    if f~=-1,fclose(f);end
    otf=outf{i};
    lst=strfind(otf,'.')-1;
    if isempty(lst),lst=length(otf);end
    otf=otf(1:lst);
    csvwrite([otf,'_user_hist.csv'],smp_hist{i});
    csvwrite([otf,'_sim_hist.csv'],sim_hist{i});
end

function out=batch_comp_encode(input_smpd,ip_smpd,inputf,ipf,input_smp_id,ip_smp_id,exp_id,bld,outf)
out=0;
s=cell(length(outf),1);
odl=s;pl=s;
load(['mm9_sn_models.mat']);
mm9_tf_beds=tf_beds;mm9_tf_dists=tf_dists;
load(['hg19_sn_models.mat']);
hg19_tf_beds=tf_beds;hg19_tf_dists=tf_dists;
clear tf_beds tf_dists;
for i=1:length(s)
    if strcmp(bld{i},'hg19'),tf_beds{i}=hg19_tf_beds;tf_dists{i}=hg19_tf_dists;
    else,tf_beds{i}=mm9_tf_beds;tf_dists{i}=mm9_tf_dists;end
end
if ~matlabpool('size'),matlabpool;end
parfor i=1:length(s)
    t=cell(3,1);
    input_sample=containers.Map;ip_sample=containers.Map;
    input_sample=input_smpd{i};ip_sample=ip_smpd{i};
    ip_data=ip_sample(ip_smp_id{i});input_data=input_sample(input_smp_id{i});
    [od,p]=find_tf_binding_odds(exp_id{i},ip_data.dens,input_data.dens,1000,tf_beds{i},tf_dists{i});
    odl{i}=od;pl{i}=p;
    t{1}=['Signal to noise ratio (SNR): ' num2str(od)];
    t{2}=['The probability of observing the given SNR or less in the ENCODE database: ' num2str(p)];
    t{3}=['A small probability indicates your data differs greatly from ENCODE datasets'];
    s{i}=t;
end
if matlabpool('size'),matlabpool close;end
for i=1:length(s)
    try, f=fopen(outf{i},'w');catch me, disp(['error opening output file ' outf{i}]),end
    fprintf(f,'IP_file,%s\n',ipf{i});fprintf(f,'Input_file,%s\n',inputf{i});
    fprintf(f,'IP_sample_id,%s\n',ip_smp_id{i});fprintf(f,'Input_sample_id,%s\n',input_smp_id{i});
    fprintf(f,'experiment_id,%s\n',exp_id{i});fprintf(f,'build,%s\n',bld{i});
    fprintf(f,'odds_ratio,%g\n',odl{i});fprintf(f,'probability,%g\n',pl{i});
    t=s{i};for j=1:length(t),fprintf(f,'%s\n',t{j});end
    if f~=-1,fclose(f);end
end

function out=batch_ip_strength(input_smpd,ip_smpd,inputf,ipf,input_smp_id,ip_smp_id,outf)
out=0;
s=cell(length(outf),1);%s{i} is a cell array of strings holding the result of the test to be stored in outf{i}
fdl=s;htl=s;kl=s;ml=s;sz_ipl=s;sz_inputl=s;pl=s;ql=s;
if ~matlabpool('size'),matlabpool;end
parfor i=1:length(s)
    input_sample=containers.Map;ip_sample=containers.Map;
    input_sample=input_smpd{i};ip_sample=ip_smpd{i};
    [tmps,fd,ht,k,m,sz_ip,sz_input,p,q]=ip_strength(input_sample(input_smp_id{i}),ip_sample(ip_smp_id{i}));
    s{i}=tmps;fdl{i}=fd;htl{i}=ht;kl{i}=k;ml{i}=m;sz_ipl{i}=sz_ip;
    sz_inputl{i}=sz_input;pl{i}=p;ql{i}=q;
end
if matlabpool('size'),matlabpool close;end
for i=1:length(s)
    try, f=fopen(outf{i},'w');catch me, disp(['error opening output file ' outf{i}]),end
    fprintf(f,'IP_file,%s\n',ipf{i});fprintf(f,'Input_file,%s\n',inputf{i});
    fprintf(f,'IP_sample_id,%s\n',ip_smp_id{i});fprintf(f,'Input_sample_id,%s\n',input_smp_id{i});
    fprintf(f,'pass,%g\n',htl{i});fd=fdl{i};
    fprintf(f,'fdr,%g\n',fd('all'));fprintf(f,'tfbs_normal_fdr,%g\n',fd('tfbs_normal'));
    fprintf(f,'histone_normal_fdr,%g\n',fd('histone_normal'));fprintf(f,'tfbs_cancer_fdr,%g\n',fd('tfbs_cancer'));
    fprintf(f,'histone_cancer_fdr,%g\n',fd('histone_cancer'));fprintf(f,'percent_genome_enriched,%g\n',(100-100*kl{i}/ml{i}));
    fprintf(f,'input_scaling_factor,%g\n',(pl{i}*sz_ipl{i})/(ql{i}*sz_inputl{i}));
    fprintf(f,'differential_percentage_enrichment,%g\n',100*(ql{i}-pl{i}));
    t=s{i};for j=1:length(t),fprintf(f,'%s\n',t{j});end
    if f~=-1,fclose(f);end
end

function [t,fd,ht,k,m,sz_ip,sz_input,p,q]=ip_strength(input_data,ip_data)
    t={};
    [p,q,ht,pval,k,m,sz_ip,sz_input,f,err]=comp_scaling_factor(ip_data.dens,input_data.dens,[0,0]);
    t{1}=['p,' num2str(p)];t{2}=['q,' num2str(q)];t{3}=['ht,' num2str(ht)];
    t{4}=['pval,' num2str(pval)];t{5}=['k,' num2str(k)];t{6}=['m,' num2str(m)];
    t{7}=['sz_ip,' num2str(sz_ip)];t{8}=['sz_input,' num2str(sz_input)];
    err_str={};err_idx=1;
    if any(err==1)
        t{length(t)+1}='';    
        t{length(t)+1}='The IP channel is extremely zero-inflated,';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='indicating a potentially insufficient depth';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='of sequencing in the IP channel.';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='The false negative rate in peak calling';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='may be high as a result.';    
        err_str{err_idx}=t{end};err_idx=err_idx+1;
    end
    if any(err==2)
        t{length(t)+1}='';
        t{length(t)+1}='The Input channel is extremely zero-inflated,';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='indicating a potentially insufficient depth';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='of sequencing in the Input channel.';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='The false positive rate in peak calling';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='may be high as a result.';    
        err_str{err_idx}=t{end};err_idx=err_idx+1;
    end
    if any(err==4)
        t{length(t)+1}='';
        t{length(t)+1}='Possible PCR amplification bias in Input,';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='more than 25% of the reads map to less';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='than 0.01% of the genome. Consider';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='de-duplicating your reads and';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='running CHANCE again.';    
        err_str{err_idx}=t{end};err_idx=err_idx+1;
    end
    if any(err==3)
        t{length(t)+1}='';
        t{length(t)+1}='The IP appears weak';  
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='The Input channel shows greater enrichment';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        t{length(t)+1}='for signal than the IP channel';
        err_str{err_idx}=t{end};err_idx=err_idx+1;
        fd=containers.Map;
        fd('all')=1;fd('tfbs_normal')=1;
        fd('histone_normal')=1;fd('tfbs_cancer')=1;fd('histone_cancer')=1;
    else
        try
        load('fdr_data.mat','qval');
        fd=containers.Map;
        kz=qval.keys;
        fdflg=1;
        for i=1:length(kz)
            %identify the fdr corresponding to the p-value for the user's data
            if pval>=max(qval(kz{i}).t),fd(kz{i})=max(qval(kz{i}).q);ht=0;
            elseif pval<=min(qval(kz{i}).t), fd(kz{i})=min(qval(kz{i}).q);ht=1;
            else
                fd(kz{i})=interp1(qval(kz{i}).t,qval(kz{i}).q,pval);
                fdflg=(fdflg && (fd(kz{i})>0.05|isnan(fd(kz{i}))|isinf(fd(kz{i}))));
            end
        end
        if fdflg,ht=0;end
        if ht==0,
            out_str={};out_idx=1;
            t{length(t)+1}='';
            t{length(t)+1}='The IP appears weak.';
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}='The percentage enrichment of IP over';
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}='Input is not statistically significant';
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=['Divergence test q-value (pFDR) is high in all samples: '];
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=sprintf('all samples\ttfbs normal\thistone normal\ttfbs cancer\thistone cancer');
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=sprintf('%f\t%f\t%f\t%f\t%f',[fd('all') ...
                                fd('tfbs_normal') fd('histone_normal') ...
                                fd('tfbs_cancer') fd('histone_cancer')]);
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            if ~isempty(err_str),out_str{out_idx}='Errors:';out_idx=out_idx+1;end
            for i=1:length(err_str),out_str{out_idx}=err_str{i};out_idx=out_idx+1;end
        else
            out_str={'IP appears successful'};out_idx=2;
            t{length(t)+1}='';
            t{length(t)+1}='Significant enrichment for signal in IP over Input';
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}='Divergence test q-value (pFDR): ';
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=sprintf('all samples\ttfbs normal\thistone normal\ttfbs cancer\thistone cancer');
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=sprintf('%f\t%f\t%f\t%f\t%f',[fd('all') ...
                                fd('tfbs_normal') fd('histone_normal') ...
                                fd('tfbs_cancer') fd('histone_cancer')]);
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=['Scaling factor: ' num2str((p*sz_ip)/(q*sz_input)) ' (scale input by this amount)'];
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=[num2str(100-100*k/m) '% of the genome is enriched for signal.'];
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            t{length(t)+1}=['Cumulative % enrichment of IP over Input: ' num2str(100*(q-p)) '%'];
            out_str{out_idx}=t{end};out_idx=out_idx+1;
            if ~isempty(err_str),out_str{out_idx}='Warnings:';out_idx=out_idx+1;end
            for i=1:length(err_str),out_str{out_idx}=err_str{i};out_idx=out_idx+1;end
        end
        catch me
            save('error.mat','me')
            alert('error computing IP strength, see error.mat')
        end
    end

function smpd=par_bin_data(fin,smp_id,bld,chr_lens,typ)
if ~matlabpool('size'),matlabpool;end
smpd={};
parfor i=1:length(fin)
    [d,nuc_freq,phred_hist,~,~]=make_density_from_file(fin{i},chr_lens{i},1000,typ{i});
    if isempty(d)
        fe=fopen('read_errors.txt','a');
        fprintf(fe,'%s\n',fin{i});
        fclose(fe);
    else
        chrs=d.keys;n=0;
        for j=1:length(chrs),n=n+sum(d(chrs{j}));end
        smp=struct('nreads',n,'genome',bld{i},'dens',d,'nuc_freq',nuc_freq,'phred',phred_hist);
        smpd{i}=containers.Map(smp_id{i},smp);
    end
end
if matlabpool('size'),matlabpool close;end

function out=disp_help()
s=sprintf('CHANCE usage:\n');
s=[s,sprintf('run_chance.sh /PATH/TO/MCR binData -b build -t file_type -s sample_id -o output_file -f file\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR binData -p parameters_file\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR IPStrength -b build -t file_type -o output_file --ipfile IP_file_name (--ipsample IP_sample_name) --inputfile input_file_name (--inputsample input_sample_name)\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR IPStrength -p parameters_file\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR multiIPNorm -p parameters_file\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR compENCODE -b build -t file_type -o output_file -e experiment_type --ipfile IP_file_name (--ipsample IP_sample_name) --inputfile input_file_name (--inputsample input_sample_name)\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR compENCODE -p parameters_file\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR spectrum -b build -t file_type (-s sample_id) -o output_file -f file_name -s sample_id\n')];
s=[s,sprintf('run_chance.sh /PATH/TO/MCR spectrum -p parameters_file\n')];
disp(s);
out=0; 