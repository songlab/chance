function [p,q,ht,pval,k,m,err] = extract_sig(s1,s2,a1,a2)
%function [p,q,ht,pval,k,m,err] = extract_sig(s1,s2,a1,a2)
%
%IN: s1 and s2 are vectors containing the alignment densities of two
%samples, the order stats are based on a sorting of s1 so if
%comparing IP to Input s1 should be IP
%    a1 is a axes handle for the Lorenz plot, a2 is a axes handle for the linearization plot  
%
%Out: 0<=p,q<=1 are percentages, which give the percentage allocation 
%     of tags to background in s1 and s2 respectively
%     ht is 1 if the null hypothesis of the binary entropy divergence test
%     (that the change in entropy between p and q is 0 to linear order) can
%     can be rejected at the %5 significance level, ht ==0 otherwise
%     pval is the pvalue associated with the above test
%     k is the index into the order stat of s1 separating signal bins from background bins
%     m is the total number of bins
%     err is a vector of error codes
%       if any(err==1), then the IP may be zero-inflated indicating
%         possible insufficiency of sequencing depth in the IP
%         implying a high false negative rate in peak detection
%       if any(err==2), then the Input may be zero-inflated
%         indicating possible insufficient sequencing depth in the
%         Input channel and implying a high false discovery rate in
%         peak detection
%       if any(err==3), then the IP Lorenz curve crosses the Input
%       curve from below indicating greater enrichment for signal
%       in the Input channel than the IP
%       if any(err==4), then there may be potential PCR amplification bias in the input channel,
%       less than .01% of the genome has more than 25% of the tags

debug=0;
err=[];
%compute the trimmed means of the order stat of ip and reordered
%input
m=length(s1);
[ss1,idx]=sort(s1);
s2r=s2(idx);
cs1=cumsum(ss1);
cs2=cumsum(s2r);
%cut the leading zeros from signal 1 (IP)
gz=max(find(ss1==0))+1;
ss1_cut=ss1(gz:end);
s2r_cut=s2r(gz:end);
cs1_cut=cumsum(ss1_cut);
cs2_cut=cumsum(s2r_cut);
%check for err 2, zero inflated input as indicated by the Input curve crossing
% the IP curve from below
if any(cs1((gz+1):floor(end/2))/cs1(end)>cs2((gz+1):floor(end/2))/cs2(end))
  err(end+1)=2;
  input_crs=gz+1;
  %find the crossing point
  while (input_crs<length(cs1))&&(cs1(input_crs)/cs1(end)>cs2(input_crs)/cs2(end)) 
      input_crs=input_crs+1;
  end
end
gz2=max(find(cs1_cut/cs1_cut(end)<cs2_cut/cs2_cut(end)));
%compute the point of maximal difference for the cut dataset
[cut_df,k]=max(abs(cs2_cut(1:gz2)/cs2_cut(end)-cs1_cut(1:gz2)/cs1_cut(end)));
%compute the percentage tag density, adding the zero bins back in
p=cs1(k+gz)/cs1(end);q=cs2(k+gz)/cs2(end);
%check for err 1, zero inflated IP as indicated by a global maximal distance
%between the two curves occuring when the IP curve is zero
[df,mx_k]=max(abs(cs1(1:gz-1)/cs1(end)-cs2(1:gz-1)/cs2(end)));
if df>0.75*cut_df|mx_k/m>0.5,err(end+1)=1;end
%compute the change in entropy and its statistical sig via div test
[ht,pval,ci,zval]=bin_ent_stat(p,q);
k=k+gz;
%check for err 3, failed IP, as indicated by the IP curve crossing the
%Input curve from below
if any(cs1(k:floor(end-end*.05))/cs1(end)>cs2(k:floor(end-end*.05))/cs2(end))
    err(end+1)=3;
    ip_crs=length(cs1)-1;
    while (ip_crs>1)&&(cs1(ip_crs)/cs1(end)>cs2(ip_crs)/cs2(end))
        ip_crs=ip_crs-1;
    end
end
%check for err 4, PCR amplification bias, detected by large "point spikes"
% in alignment density
k_pcr=m-floor(0.01*m);y_pcr=cs2(k_pcr)/cs2(end);
if y_pcr<.75 %set the error flag and rerun analysis discarding last 0.01% of the order stat
    err(end+1)=4;
    ss1_tcut=ss1(1:k_pcr);
    s2r_tcut=s2r(1:k_pcr);
    cs1_tcut=cumsum(ss1_tcut);
    cs2_tcut=cumsum(s2r_tcut);
    ss1_cut=ss1(gz:k_pcr);
    s2r_cut=s2r(gz:k_pcr);
    cs1_cut=cumsum(ss1_cut);
    cs2_cut=cumsum(s2r_cut);
    [~,kt]=max(abs(cs2_cut/cs2_cut(end)-cs1_cut/cs1_cut(end)));
    kt=kt+gz;
    p=cs1_tcut(kt+gz)/cs1_tcut(end);q=cs2_tcut(kt+gz)/cs2_tcut(end);
    [ht,pval,ci,zval]=bin_ent_stat(p,q);
end
if debug, keyboard();end

%lorenz plot
if ~isempty(a1)
    m=length(cs1);
    plot(a1,linspace(0,1,m),cs1/cs1(end),'b','LineWidth',3)
    hold(a1,'on')
    plot(a1,linspace(0,1,m),cs2/cs2(end),'r','LineWidth',3)
    title(a1,'Cumulative percentage enrichment in each channel','FontName','Arial','FontSize',20)
    xlabel(a1,'Percentage of bins','FontName','Arial','FontSize',20)
    ylabel(a1,'Percentage of tags','FontName','Arial','FontSize',20)
    set(a1,'FontName','Arial','FontSize',20)
    lgnd_str={'IP','Input'};
    if any(err==3)
        plot(a1,[ip_crs/m ip_crs/m],[0 1],'g--','LineWidth',3)
        lgnd_str{end+1}=sprintf('Possible failed IP');
    else
        plot(a1,[k/m k/m],[0 1],'g--','LineWidth',3);
        lgnd_str{end+1}=sprintf('Signal/Background cutoff');
    end  
    if any(err==1)
        plot(a1,[mx_k/m mx_k/m],[0 1],'k--','LineWidth',3);
        lgnd_str{end+1}=sprintf('Insufficient sequencing\ndepth in the IP channel');
    end
    if any(err==2)
        plot(a1,[input_crs/m input_crs/m],[0 1],'m--','LineWidth',3)    
        lgnd_str{end+1}=sprintf('Insufficient sequencing\ndepth in the Input channel');
    end
    if any(err==4)
        plot(a1,[0 1],[y_pcr y_pcr],'c--','LineWidth',3)
        lgnd_str{end+1}=sprintf('Possible amplification bias');
    end
    legend(a1,lgnd_str,'Location','NorthWest')
    hold(a1,'off')
end
%linearization plot
if ~isempty(a2)
    y=(cs1./cs2)*(cs2(end)/cs1(end));y=y(1:1000:end);
    t1=linspace(0,1,length(y))';
    q1=quantile(y,[.30 .70]);idx1=find(y>0&y>=min(q1)&y<=max(q1));
    plot(a2,t1(idx1),y(idx1),'o')
    hold(a2,'on')
    ck=[max(idx1)+1:length(y)];
    scatter(a2,t1(ck),y(ck),'ro')
    lgnd=legend('30th-70th quantiles', 'Ratio of trimmed means','Location','North');
    set(lgnd,'FontName','Arial','FontSize',20)
    axis tight;
    %if ~any(err),plot(a2,[k/m k/m],[0 1],'g');end
    title('Linearization plot','FontName','Arial','FontSize',20)
    xlabel(a1,'Percentage of bins','FontName','Arial','FontSize',20)
    set(a1,'FontName','Arial','FontSize',20)

end