function varargout = chip_qc_main(varargin)
% CHIP_QC_MAIN MATLAB code for chip_qc_main.fig
%      CHIP_QC_MAIN, by itself, creates a new CHIP_QC_MAIN or raises the existing
%      singleton*.
%
%      H = CHIP_QC_MAIN returns the handle to a new CHIP_QC_MAIN or the handle to
%      the existing singleton*.
%
%      CHIP_QC_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHIP_QC_MAIN.M with the given input arguments.
%
%      CHIP_QC_MAIN('Property','Value',...) creates a new CHIP_QC_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chip_qc_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chip_qc_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chip_qc_main

% Last Modified by GUIDE v2.5 21-May-2012 10:23:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chip_qc_main_OpeningFcn, ...
                   'gui_OutputFcn',  @chip_qc_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before chip_qc_main is made visible.
function chip_qc_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chip_qc_main (see VARARGIN)

splash('chance.png');


% Choose default command line output for chip_qc_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chip_qc_main wait for user response (see UIRESUME)
% uiwait(handles.root_window);

% --- Outputs from this function are returned to the command line.
function varargout = chip_qc_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes during object creation, after setting all properties.
function main_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',{'Welcome to CHANCE!'});

% --- Executes on button press in load_sample.
function load_sample_Callback(hObject, eventdata, handles)
% hObject    handle to load_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
%sample_data is a Map between sample id strings and data objects
%each data object has a density Map dens, a phred histogram phred,
%a nucleotide frequency histogram nuc_freq, and the number of reads nreads
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),sample_data=containers.Map;end
if ~isfield(main_data,'last_dir')
    [fname,pname,filteridx]=uigetfile('*.*','Choose a file.');
    if fname==0,return;end
else
    [fname,pname,filteridx]=uigetfile([main_data.last_dir '*.*'],'Choose a file.');
    if fname==0,return;end
end
ot=choose_file_type('title','Choose file type and genome','Select file type and build.');
if strcmp(ot.genome,'other')
    if ~(strcmp(ot.file_type,'bam')|strcmp(ot.file_type,'sam'))
        alert('String','For organisms other that human or mouse the file type must be SAM or BAM')
        return;
    end
    alert('String','Note: any CHANCE functions or statistics dependent on ENCODE data are only available for human or mouse.');
end
main_data.last_dir=pname;
if ~isfield(main_data,'chrom_lens')||~strcmp(main_data.genome,ot.genome)    
    main_data.genome=ot.genome;
    if ~strcmp(ot.file_type,'bam')
        load([ot.genome 'lengths.mat']);
        main_data.chrom_lens=chr_lens;
    else
        main_data.chrom_lens=[];
    end
end
set(handles.root_window,'UserData',main_data);
t=get(handles.main_output,'String');
t{length(t)+1}=['Reading ' fname '... this may take a few minutes.'];
set(handles.main_output,'String',t);
pause(0.5)
try
    [d,nuc_freq,phred,chr_lens,cncl]=make_density_from_file([pname fname],main_data.chrom_lens,1000,ot.file_type);
catch me
    alert('title','File read error','string','Make sure file matches selected type.');
    t=get(handles.main_output,'String');
    t{length(t)+1}=['Error reading ' fname ' make sure file matches selected file type.'];
    set(handles.main_output,'String',t);
    return
end
if cncl
    alert('title','File processing canceled','string',[fname ' processing canceled.']);
    t=get(handles.main_output,'String');
    t{length(t)+1}=[fname ' processing canceled.'];
    set(handles.main_output,'String',t);
else
    main_data.chrom_lens=chr_lens;
    set(handles.root_window,'UserData',main_data);
    try
        smp_name=set_sample_id('title','Enter sample ID:','string',sprintf(['Enter a name for the sample\n(' fname ')']));
    catch me
        smp_name='new_sample';
    end
    chrs=d.keys;n=0;
    for i=1:length(chrs),n=n+sum(d(chrs{i}));end
    smp.nreads=n;smp.genome=ot.genome;
    smp.dens=d;smp.nuc_freq=nuc_freq;smp.phred=phred;
    smp_id=[smp_name ' - # reads: ' num2str(n) ' - genome: ' smp.genome ' - ' fname];
    sample_data(smp_id)=smp;
    smp_strs=cellstr(get(handles.sample_list,'String'));
    if strcmp(smp_strs{1},'No sample loaded...'),smp_strs={};end
    smp_strs{length(smp_strs)+1}=smp_id;
    set(handles.sample_list,'String',smp_strs);
    set(handles.load_sample,'UserData',sample_data);
    t=get(handles.main_output,'String');
    t{length(t)+1}=[fname ' processed successfully, ' num2str(n) ' tags processed.'];
    set(handles.main_output,'String',t);
end 

% --- Executes on button press in nuc_freq.
function nuc_freq_Callback(hObject, eventdata, handles)
% hObject    handle to nuc_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
smp_id=contents{get(handles.sample_list,'Value')};
%smp_id=sscanf(smp_id,'%[^-]');smp_id=smp_id(1:end-1);
dt=sample_data(smp_id);
f1=figure;
a1=gca;
hold(a1,'on')
nf=dt.nuc_freq;
plot(a1,nf.A,'b','LineWidth',3)
plot(a1,nf.C,'g','LineWidth',3)
plot(a1,nf.G,'r','LineWidth',3)
plot(a1,nf.T,'y','LineWidth',3)
set(f1,'color','w')
set(a1,'FontName','Arial','FontSize',20)
stp=strfind(smp_id,' - # reads')-1;if isempty(stp),stp=min(length(smp_id),15);end
title(a1,[smp_id(1:stp) ' nucleotide frequencies'],'FontName','Arial','FontSize',20)
xlabel(a1,'Base position','FontName','Arial','FontSize',20)
ylabel(a1,'Frequency','FontName','Arial','FontSize',20)
l=legend(a1,'A','C','G','T');
set(l,'FontName','Arial','FontSize',20);

f2=figure;
a2=gca;
plot(a2,nf.N,'b','LineWidth',3)
set(f2,'color','w')
set(a2,'FontName','Arial','FontSize',20)
stp=strfind(smp_id,' - # reads')-1;if isempty(stp),stp=min(length(smp_id),15);end
title(a2,['Frequency of uncallable bases in ' smp_id(1:stp)],'FontName','Arial','FontSize',20)
xlabel(a2,'Base position','FontName','Arial','FontSize',20)
ylabel(a2,'Frequency','FontName','Arial','FontSize',20)

% --- Executes on button press in phred_score.
function phred_score_Callback(hObject, eventdata, handles)
% hObject    handle to phred_score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
smp_id=contents{get(handles.sample_list,'Value')};
%smp_id=sscanf(smp_id,'%[^-]');smp_id=smp_id(1:end-1);
dt=sample_data(smp_id);
ps=choose_phred();%get the phred score offset based on sequencing convention
cs=flipud(cumsum(dt.phred));
for i=1:size(cs,2),cs(:,i)=cs(:,i)/cs(1,i);end%comp cumulative phred score densities
k=min(find(sum(cs')==0))+10;
p=max(find(size(cs,2)==sum(cs')))-10;
cs=1-cs(p:k,:);
f1=figure;a1=gca;
heatmap(cs,1:size(cs,2),(size(cs,1)+k-p-ps):-1:(1+k-p-ps));
set(f1,'color','w')
set(a1,'FontName','Arial','FontSize',20)
stp=strfind(smp_id,' - # reads')-1;if isempty(stp),stp=min(length(smp_id),15);end
title(a1,[smp_id(1:stp) ' Phred quality score cumulative density'],'FontName','Arial','FontSize',20)
colorbar('peer',a1,'FontName','Arial')
xlabel(a1,'Base position','FontName','Arial','FontSize',20)
ylabel(a1,'Quality score','FontName','Arial','FontSize',20)

% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile('*.txt','Save data as...','chipseq_qc_log.txt');
else
    [fname pname]=uiputfile([main_data.last_dir '*.txt'],'Save data as...','chipseq_qc_log.mat');
end
main_data.last_dir=pname;
set(handles.root_window,'UserData',main_data);
f=fopen([pname fname],'w');
s=get(handles.main_output,'String');
fprintf(f,'%s\n',s{:});
fclose(f);
alert('title','Wrote session log','string',['Wrote log file ' fname]);
t=get(handles.main_output,'String');
t{length(t)+1}=['Wrote log file ' fname];
set(handles.main_output,'String',t);


% --- Executes on button press in extract_signal.
function extract_signal_Callback(hObject, eventdata, handles)
% hObject    handle to extract_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


plot_on=[1,1];
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
ip_id=contents{get(handles.sample_list,'Value')};
%ip_id=sscanf(ip_id,'%[^-]');ip_id=ip_id(1:end-1);
ip_data=sample_data(ip_id);
try
    input_id=choose_sample('title','Chose Input:','string','Choose an Input (control) sample...','smp_lst',contents);
catch me
    return
end
%input_id=sscanf(input_id,'%[^-]');input_id=input_id(1:end-1);
input_data=sample_data(input_id);
if ~strcmp(ip_data.genome,input_data.genome)
    alert('title','Genome mismatch!','String','The genomes of the IP and Input don''t match')
    return
end
[p,q,ht,pval,k,m,sz_ip,sz_input,f,err]=comp_scaling_factor(ip_data.dens,input_data.dens,plot_on);
t=get(handles.main_output,'String');
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
        f1=figure;
        set(f1,'color','w');
        ax=gca;set(ax,'Visible','off');
        text(.5,.5,out_str,'FontName','Arial','HorizontalAlignment','center','FontSize',20)
        title(ax,'Summary','FontName','Arial','FontSize',20);
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
        f1=figure;
        set(f1,'color','w');
        ax=subplot(2,2,1);
        pie(ax,[1-q+p q-p],{sprintf(['%%' num2str(100*(1-q+p)) '\nbackground reads']),sprintf(['%%' num2str(100*(q-p)) '\nsignal reads'])})
        tx=findobj(ax,'Type','text');
        set(tx,'FontName','Arial','FontSize',18)
        title(ax,'IP strength','FontName','Arial','FontSize',18);
        ax=subplot(2,2,2);
        pie(ax,[k/m (1-k/m)],{sprintf(['%%' num2str(100*k/m) ' is\ndevoid of signal']),sprintf(['%%' num2str(100-100*k/m) ' of genome is\nenriched for signal'])})
        tx=findobj(ax,'Type','text');
        set(tx,'FontName','Arial','FontSize',18)
        ax=subplot(2,2,[3 4]);set(ax,'Visible','off');
        text(.5,.5,out_str,'FontName','Arial','HorizontalAlignment','center','FontSize',18)
        title(ax,'Summary','FontName','Arial','FontSize',18);
        colormap(f1,'cool')
    end
    catch me
        save('error.mat','me')
        alert('problems','problems')
    end
end
set(handles.main_output,'String',t);



% --- Executes on button press in comp_spectrum.
function comp_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to comp_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
smp_id=contents{get(handles.sample_list,'Value')};
%sample_id=sscanf(sample_id,'%[^-]');sample_id=sample_id(1:end-1);
sample_data=sample_data(smp_id);
Smpl=[];dens=sample_data.dens;chrs=dens.keys;
for i=1:length(chrs),Smpl=[Smpl;dens(chrs{i})];end
[c,l]=wavedec(Smpl,15,'haar');
[ea,ed_inp]=wenergy(c,l);
t=get(handles.main_output,'String');
t{length(t)+1}='';
t{length(t)+1}=['Energy in approximation coifficients: ' num2str(ea)];
set(handles.main_output,'String',t);
d=fitdist(max(Smpl,1),'gamma');
sim_data=poissrnd(d.random(length(Smpl),1));
[c,l]=wavedec(sim_data(find(~isnan(sim_data))),15,'haar');
[ea,ed_sim]=wenergy(c,l);
f1=figure;a1=gca;
set(f1,'color','w');
bar(a1,[ed_inp'/sum(ed_inp),ed_sim'/sum(ed_sim)])
stp=strfind(smp_id,' - # reads')-1;if isempty(stp),stp=min(length(smp_id),15);end
legend(a1,smp_id(1:stp),'Poisson simulation')
set(a1,'FontName','Arial','FontSize',20)
xlabel(a1,'Characteristic length scale in kbp','FontName','Arial','FontSize',20);
for i=1:15,tl{i}=num2str(2^(i-1));end
set(a1,'XTickLabel',tl)
rotateXLabels(a1,45)
ylabel(a1,'% of spectral energy (% of variance)','FontName','Arial','FontSize',20);
title(a1,['Distribution of spectral energy in ' smp_id(1:stp)],'FontName','Arial','FontSize',20);



% --- Executes on button press in sig_noise_ratio.
function sig_noise_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to sig_noise_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
ip_id=contents{get(handles.sample_list,'Value')};
%ip_id=sscanf(ip_id,'%[^-]');ip_id=ip_id(1:end-1);
ip_data=sample_data(ip_id);
if strcmpi(ip_data.genome,'hg18')|strcmpi(ip_data.genome,'tair10')
    alert('title','Build not supported','string',[ip_data.genome ' is not supported']);
    return
end
try
    input_id=choose_sample('title','Chose Input:','string','Choose an Input (control) sample...','smp_lst',contents);
catch me
    return
end
%input_id=sscanf(input_id,'%[^-]');input_id=input_id(1:end-1);
input_data=sample_data(input_id);
load([ip_data.genome '_sn_models.mat']);
try
    tfname=choose_tfname('tflist',intersect(tf_beds.keys,tf_dists.keys));
catch me
    return
end
[od,p,h,a1]=find_tf_binding_odds(tfname,ip_data.dens,input_data.dens,1000,tf_beds,tf_dists);
set(h,'Visible','on');
set(h,'color','w');
set(a1,'FontName','Arial','FontSize',20);
xlabel(a1,'log2 odds','FontName','Arial','FontSize',20);
set(a1,'YTick',[]);
title(a1,'log2 IP/Input odds ratio a random tag lies in a consensus peak','FontName','Arial','FontSize',20);
legend(a1,'normal model','ENCODE data','your experiment');
t=get(handles.main_output,'String');
t{length(t)+1}='';
t{length(t)+1}=['Signal to noise ratio (SNR): ' num2str(od)];
t{length(t)+1}=['The probability of observing the given SNR or less in the ENCODE database: ' num2str(p)];
t{length(t)+1}=['A small probability indicates your data differs greatly from ENCODE datasets'];
set(handles.main_output,'String',t);



% --- Executes on button press in export_binary.
function export_binary_Callback(hObject, eventdata, handles)
% hObject    handle to export_binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.load_sample,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile('*.mat','Save data as...','Untitled.mat');
else
    [fname pname]=uiputfile([main_data.last_dir '*.mat'],'Save data as...','Untitled.mat');
end
main_data.last_dir=pname;
set(handles.root_window,'UserData',main_data);
save([pname fname],'sample_data');
alert('title','Working samples saved','string',['Wrote matlab binary file ' fname]);
t=get(handles.main_output,'String');
t{length(t)+1}=['Wrote matlab binary file ' fname];
set(handles.main_output,'String',t);


% --- Executes on button press in validate_regions.
function validate_regions_Callback(hObject, eventdata, handles)
% hObject    handle to validate_regions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
ip_id=contents{get(handles.sample_list,'Value')};
%ip_id=sscanf(ip_id,'%[^-]');ip_id=ip_id(1:end-1);
ip_data=sample_data(ip_id);
try
    input_id=choose_sample('title','Chose Input:','string','Choose an Input (control) sample...','smp_lst',contents);
catch me
    return
end
%input_id=sscanf(input_id,'%[^-]');input_id=input_id(1:end-1);
input_data=sample_data(input_id);
[p,q,ht,pval,k,m,sz_ip,sz_input,f,err]=comp_scaling_factor(ip_data.dens,input_data.dens,[0 0]);
load('div_fdr.mat');
if pval>=max(pvs),fd=max(fdrs);ht=0;
elseif pval<=min(pvs)&~any(err==3), fd=min(fdrs);ht=1;
else
    fd=interp1(pvs,fdrs,pval);
    if fd>0.05|isnan(fd)|isinf(fd)|any(err==3),ht=0;end
end
if ht==0
    t=get(handles.main_output,'String');
    t{length(t)+1}=['Warning! the IP may have failed. An accurate scaling factor'];
    t{length(t)+1}=['can not be computed. Using sequencing depth scaling.'];
    set(handles.main_output,'String',t);
    scaling_factor=sz_ip/sz_input;
else
    scaling_factor=(p*sz_ip)/(q*sz_input);
end
[fname,pname]=uigetfile([main_data.last_dir '*.*'],'Choose bed file to validate');
if fname==0,return;end
%read in bed file
f=fopen([pname fname]);
D=textscan(f,'%s%n%n%s%*[^\n]');
t=get(handles.main_output,'String');
t{length(t)+1}=[fname ' processed successfully.'];
set(handles.main_output,'String',t);
%D=={chrom_id,start,stop,region_id}
ip_tc=zeros(length(D{1}),1);
input_tc=zeros(length(D{1}),1);
pvals=ones(length(D{1}),1);
fc=pvals;
t=get(handles.main_output,'String');
t{length(t)+1}='';
t{length(t)+1}=['Region ID, IP tag count, Input tag count, Fold change, Poisson p-value'];
set(handles.main_output,'String',t);
for i=1:length(ip_tc)
   [ip_tc(i),~,~]=sum_on_intervs(ip_data.dens(D{1}{i}),[D{2}(i) D{3}(i)],1000);
   [input_tc(i),~,~]=sum_on_intervs(scaling_factor*input_data.dens(D{1}{i}),[D{2}(i) D{3}(i)],1000);
   fc(i)=ip_tc(i)/input_tc(i);
   pvals(i)=max(1-poisscdf(ip_tc(i),input_tc(i)),1e-300);
   t=get(handles.main_output,'String');
   t{length(t)+1}=[D{4}{i} ', ' num2str(ip_tc(i)) ', ' num2str(input_tc(i)) ', ' num2str(fc(i)) ', ' num2str(pvals(i))];
   set(handles.main_output,'String',t);
end
col_labels={};
for i=1:length(ip_tc)
    col_labels{i}=D{4}{i};% ' (FC:' num2str(fc(i)) ',P:' num2str(pvals(i)) ')'];
end
f1=figure;a1=gca;
set(f1,'color','w');
bar(a1,[ip_tc,input_tc])
set(a1,'xticklabel',col_labels)
rotateXLabels(a1,45)
set(a1,'FontName','Arial','FontSize',16)
legend(a1,'IP tag count','Input tag count (normalized)')
f2=figure;a2=gca;
set(f2,'color','w');
out_str={'Validation results:',''};out_idx=3;
for i=1:length(ip_tc)
    out_str{out_idx}=[D{4}{i} ', IP/Input fold change:' num2str(fc(i)) ', p-value:' num2str(pvals(i))];
    out_idx=out_idx+1;
end
set(a2,'Visible','off');
text(.5,.5,out_str,'FontName','Arial','HorizontalAlignment','center','FontSize',18)
title(a2,'Validation results','FontName','Arial','FontSize',18);



% --- Executes on button press in restore_session.
function restore_session_Callback(hObject, eventdata, handles)
% hObject    handle to restore_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uigetfile('*.mat','Restore session from file...');
else
    [fname pname]=uigetfile([main_data.last_dir '*.mat'],'Restore session from file...');
end
if ~isstr(fname),return;end
load([pname fname],'sample_data');
t=get(handles.main_output,'String');
t{length(t)+1}=[fname ' loaded'];
set(handles.main_output,'String',t);
if ~exist('sample_data')
    alert('title','No samples found!','string',['No samples found in ' fname])
    t=get(handles.main_output,'String');
    t{length(t)+1}=['No samples found in ' fname]; 
    set(handles.main_output,'String',t);
    main_data.last_dir=pname;
else
    kz=sample_data.keys;
    val=sample_data(kz{1});
    t=get(handles.main_output,'String');
    t{length(t)+1}=[num2str(length(kz)) ' samples loaded']; 
    set(handles.main_output,'String',t);
    main_data.last_dir=pname;
    old_sample_data=get(handles.load_sample,'UserData');
    if isempty(old_sample_data),old_sample_data=containers.Map;end
    contents = cellstr(get(handles.sample_list,'String'));
    if length(contents)==1
        if strcmp('No sample loaded...',contents{1}),contents={};end 
    end
    for i=1:length(kz)
        old_sample_data(kz{i})=sample_data(kz{i});
        contents{length(contents)+1}=kz{i};
    end
    set(handles.load_sample,'UserData',old_sample_data);
    set(handles.sample_list,'String',contents);
end
set(handles.root_window,'UserData',main_data);


% --- Executes on selection change in sample_list.
function sample_list_Callback(hObject, eventdata, handles)
% hObject    handle to sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sample_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sample_list

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(hObject,'String'));
smp_id=contents{get(hObject,'Value')};
%sample_id=sscanf(sample_id,'%[^-]');sample_id=sample_id(1:end-1);
smp=sample_data(smp_id);
if ~isfield(main_data,'chrom_lens')||~strcmp(main_data.genome,smp.genome)
    main_data.genome=smp.genome;
    load([smp.genome 'lengths.mat']); %loads chr_lens
    main_data.chrom_lens=chr_lens;
end
set(handles.root_window,'UserData',main_data);


% --- Executes during object creation, after setting all properties.
function sample_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in delete_sample.
function delete_sample_Callback(hObject, eventdata, handles)
% hObject    handle to delete_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
rmidx=get(handles.sample_list,'Value');
to_rm=contents{rmidx};
%to_rm=sscanf(to_rm,'%[^-]');to_rm=to_rm(1:end-1);
if sample_data.isKey(to_rm),sample_data.remove(to_rm);end
set(handles.load_sample,'UserData',sample_data);
k=1;
if length(contents)==1,set(handles.sample_list,'String','No sample loaded...');return,end
for i=1:length(contents)
    if ~strcmp(contents{i},contents{rmidx}), cnts2{k}=contents{i};k=k+1;end
end
if ~exist('cnts2')
    cnts2={'No sample loaded...'};
    set(handles.sample_list,'String',cnts2);
    set(handles.sample_list,'Value',1);
else
    set(handles.sample_list,'String',cnts2);
    set(handles.sample_list,'Value',1);
    smp_id=contents{get(handles.sample_list,'Value')};
    %smp_id=sscanf(smp_id,'%[^-]');smp_id=smp_id(1:end-1);
    smp=sample_data(smp_id);
    if ~isfield(main_data,'chrom_lens')||~strcmp(main_data.genome,smp.genome)    
        main_data.genome=smp.genome;
        load([smp.genome 'lengths.mat']);
        main_data.chrom_lens=chr_lens;
    end
end


function multi_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to multi_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%identify the samples to work with
sample_data=get(handles.load_sample,'UserData');
if isempty(sample_data),return;end
contents = cellstr(get(handles.sample_list,'String'));
num_samples=get_num_samples('title','Set number of samples...','String','How many IP samples will you compare?');
ip_samples={};rS=[];
for i=1:num_samples
    try
        ip_samples_id{i}=choose_sample('title','Chose an IP sample...','string',['Choose IP sample ' num2str(i)],'smp_lst',contents);
    catch me
        return
    end
    %tmp_id=sscanf(ip_samples_id{i},'%[^-]');ip_samples_id{i}=tmp_id(1:end-1);
    ip_samples{i}=sample_data(ip_samples_id{i});
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
f1=figure;a1=gca;f2=figure;set(f1,'color','w');set(f2,'color','w');
a2=gca;set(a2,'Visible','off','FontName','Arial','FontSize',20);
plot(a1,linspace(0,1,m),cs1/cs1(end),'g-.','LineWidth',3)
cc=colormap('HSV');
hold(a1,'on');
t=get(handles.main_output,'String');
t{length(t)+1}='';
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
    plot(a1,linspace(0,1,m),CS2(:,i)/CS2(end,i),'LineWidth',3,'Color',cc(floor(i*size(cc,1)/num_samples),:))
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
set(handles.main_output,'String',t);
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
%fh1=figure;ah1=gca;
%eatmap(100*diff_enrich,ip_samples_id,ip_samples_id,'%2.1f%%','FontSize',20,'TickAngle',45,'ShowAllTicks',true,'Colormap','cool');
%set(fh1,'color','w');
%set(ah1,'FontName','Arial','FontSize',20);
%title(ah1,'Cumulative differential enrichment','FontName','Arial','FontSize',20);
%colorbar('peer',ah1,'FontName','Arial','FontSize',20);
fh2=figure;ah2=gca;
for i=1:length(ip_samples_id)
    tmp_str=ip_samples_id{i};
    stp=strfind(tmp_str,' - # reads')-1;
    if isempty(stp),stp=min(length(tmp_str),10);end
    ssmp_id{i}=tmp_str(1:stp);
end
heatmap(100*diff_genom,ssmp_id,ssmp_id,'%2.1f%%','FontSize',30,'TickAngle',45,'ShowAllTicks',true,'Colormap','cool','TickFontSize',30);
set(fh2,'color','w');
set(ah2,'FontName','Arial','FontSize',30);
title(ah2,'Percent of genome differentially enriched','FontName','Arial','FontSize',30);
colorbar('peer',ah2,'FontName','Arial','FontSize',30);
title(a1,'Cumulative percentage enrichment in each channel','FontName','Arial','FontSize',20)
xlabel(a1,'Percentage of bins','FontName','Arial','FontSize',20)
ylabel(a1,'Percentage of tags','FontName','Arial','FontSize',20)
lbls={'consensus'};lbl_idx=2;
for i=1:length(ip_samples_id),lbls{lbl_idx}=ssmp_id{i};lbl_idx=lbl_idx+1;end
legend(a1,lbls, 'Location','NorthWest')
set(a1,'FontName','Arial','FontSize',20)
axes(a2);text(0.5,0.5,txt_out,'FontName','Arial','HorizontalAlignment','center','FontSize',18)
