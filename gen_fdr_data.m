function [pnull,palt,cinull,cialt,zvalnull,zvalalt]=gen_fdr_data()
addpath /home/aaron/research/fun_genom/matlab_code/chance
pref='/home/aaron/research/fun_genom/data/uchicago_tfbs/';
d=dir([pref '*.mat']);
pnull=[];palt=[];cinull=[];cialt=[];zvalnull=[];zvalalt=[];
IPZ={};INPUT1Z={};INPUT2Z={};
for i=1:length(d)
    disp(d(i).name)
    load([pref d(i).name]);
    input1_struc=sample_data('input1');
    input2_struc=sample_data('input2');
    ip_struc=sample_data('ip');
    kz=input1_struc.dens.keys;input1=[];
    for j=1:length(kz),input1=[input1;input1_struc.dens(kz{j})];end
    kz=input2_struc.dens.keys;input2=[];
    for j=1:length(kz),input2=[input2;input2_struc.dens(kz{j})];end
    kz=ip_struc.dens.keys;ip=[];
    for j=1:length(kz),ip=[ip;ip_struc.dens(kz{j})];end
    IP=ones(size(ip,1),10);INPUT1=ones(size(input1,1),10);INPUT2=ones(size(input2,1),10);
    for j=1:length(ip),IP(j,:)=ip(j);end
    for j=1:length(input1),INPUT1(j,:)=input1(j);end
    for j=1:length(input2),INPUT2(j,:)=input2(j);end
    IP=[poissrnd(IP),ip];INPUT1=[poissrnd(INPUT1),input1];INPUT2=[poissrnd(INPUT2),input2];
    IPZ{i}=IP;INPUT1Z{i}=INPUT1;INPUT2Z{i}=INPUT2;
end
parfor i=1:length(d)
    cinull_loc=[];cialt_loc=[];
    zvalnull_loc=[];zvalalt_loc=[];
    pnull_loc=[];palt_loc=[];rn=1;ra=1;
    INPUT1=INPUT1Z{i};INPUT2=INPUT2Z{i};IP=IPZ{i};
    for j=1:11
        for k=j:11
           %compare input-input
           [p,q,~,pnull_loc(rn),~,~,~]=extract_sig(INPUT1(:,j),INPUT2(:,k),[],[]);rn=rn+1;
           [~,~,cinull_loc(:,rn),zvalnull_loc(rn)]=bin_ent_stat(p,q);
           [p,q,~,pnull_loc(rn),~,~,~]=extract_sig(INPUT2(:,k),INPUT1(:,j),[],[]);rn=rn+1;
           [~,~,cinull_loc(:,rn),zvalnull_loc(rn)]=bin_ent_stat(p,q);
           %compare ip-input
           [p,q,~,palt_loc(ra),~,~,err]=extract_sig(IP(:,j),INPUT2(:,k),[],[]);
           if any(err==3)
               palt_loc(ra)=1e-30;
               [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,p);
           else
               [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,q);
           end
           ra=ra+1;
           [p,q,~,palt_loc(ra),~,~,err]=extract_sig(IP(:,j),INPUT1(:,k),[],[]);
           if any(err==3)
               palt_loc(ra)=1e-30;
               [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,p);
           else
               [~,~,cialt_loc(:,ra),zvalalt_loc(ra)]=bin_ent_stat(p,q);
           end
           ra=ra+1;
        end
    end
    cinull=[cinull cinull_loc];zvalnull=[zvalnull zvalnull_loc];
    cialt=[cialt cialt_loc];zvalalt=[zvalalt zvalalt_loc];
    pnull=[pnull pnull_loc];palt=[palt palt_loc];
end
save uchicago_tf_cancer_div_fdr_estimates.mat