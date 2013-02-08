function out=chance_com(subr,varargin)
%function chance_com(subr,varargin)
%
%wrapper function for chance 

cmds={'binData','IPStrength','multiIPNorm','compENCODE','nucFreq','phred','spectrum'};
if ~ismember(subr,cmds),disp_help();return;end
if strcmp(subr,'binData')
    options = containers.Map({'-b','-t','-s','-o','-f'},{[],[],[],pwd,[]});
end
optionNames = options.keys;
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2,disp_help();return;end
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1});
   if any(strmatch(inpName,optionNames))
      options(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
      disp_help();
      return;
   end
end
if strcmp(subr,'binData')
    if strcmp(options('-b')
    [d,nuc_freq,phred_hist,cncl]=make_density_from_file(options('-f'),chr_lens,bin,type);
end

function out=disp_help()
s=sprintf('CHANCE usage:\n');
s=[s,sprintf('chance binData -b build -t file_type -s sample_id (-o output_directory) -f file\n')];
s=[s,sprintf('chance IPStrength (-o output_directory) --IPFile IP_file_name --IPSample IP_sample_name --InputFile input_file_name --InputSample input_sample_name\n')];
s=[s,sprintf('chance multiIPNorm -p parameters_file\n')];
s=[s,sprintf('chance compENCODE (-o output_directory) -e experiment_type --IPFile IP_file_name --IPSample IP_sample_name --InputFile input_file_name --InputSample input_sample_name\n')];
s=[s,sprintf('chance nucFreq (-o output_directory) -f file_name -s sample_id\n')];
s=[s,sprintf('chance phred (-o output_directory) -f file_name -s sample_id\n')];
s=[s,sprintf('chance spectrum (-o output_directory) -f file_name -s sample_id\n')];
disp(s);
out=0;