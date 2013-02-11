function [od,p,h,a1]=find_tf_binding_odds(tfname,ip_d,input_d,bin,tf_beds,tf_dists)
%function [od,p,h,a1]=find_tf_binding_odds(tfname,ip_d,input_d,bin,tf_beds,tf_dists)
%
%IN: fname is a string holding the name of the transcription factor to 
%      compute the odds ratio of
%    ip_d and input_d are Maps from chromosomes to alignment densities
%    bin is the number of bp which was used to bin alignments
%    tf_beds is a Map from transcription factor id strings to tf peak Maps
%    tf_dists is a Map from transcription factor id strings to structures
%      holding mean and std for a normal model of log2 odds ratios for
%      known binding sites of the given tf as well as the observed odds ratios
%      used to construct the model
%
%OUT: od is the odds ratio of read count found in known binding sites for
%       the given tf compared to input
%     h is a handle to a plot of the log2 odds ratio against the model
%     a1 is a handle to the axes the plot is into
%     p is the probability of the log odds of the given experiment or smaller
%       under the model determined by known experiments

%check that the TF is found in the database of known experiments
tf_dists_keys=tf_dists.keys;
tf_beds_keys=tf_beds.keys;
didx=find(strcmpi(tfname,tf_dists_keys));
bidx=find(strcmpi(tfname,tf_beds_keys));
if isempty(didx)|isempty(bidx)
    disp('transcription factor not found in database')
    od=1;p=1;h=0;
    return
end
[p11,~,ts]=sum_density_over_bed(ip_d,tf_beds(tf_beds_keys{bidx}));
p12=ts-p11;
[q11,~,ts]=sum_density_over_bed(input_d,tf_beds(tf_beds_keys{bidx}));
q12=ts-q11;
od=(p11*q12)/(q11*p12);
tfd=tf_dists(tf_dists_keys{didx});
p=normcdf(log2(od),tfd.mu,tfd.sigma);
h=figure('Visible','off');a1=axes;
t=linspace(tfd.mu-10*tfd.sigma,tfd.mu+10*tfd.sigma,100);
plot(a1,t,normpdf(t,tfd.mu,tfd.sigma),'LineWidth',3);
hold(a1,'on');
plot(a1,log2(tfd.ors),zeros(size(tfd.ors)),'bx','LineWidth',3);
plot(a1,log2(od),0,'r*','LineWidth',3);
