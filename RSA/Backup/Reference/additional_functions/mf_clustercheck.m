% mf_clusterchecker: function that returns which clusters in a fieldtrip
% stat variable are significant
% alpha= p-value threshold
% stat= fieldtrip stat result

function [pos_check, pos_sig, neg_check, neg_sig,sig_text]=mf_clustercheck(stat,alpha)


neg_check=0;
neg_sig=[];
sig_part{1}={'no negcluster'};
if isfield(stat, 'negclusters')
    if ~isempty(stat.negclusters)
      all_ps=[stat.negclusters.prob]
      neg_sig=find(all_ps<=alpha);
      if any(neg_sig)
    neg_check=1;
    tmp_text=repmat({'negcluster:p='},numel(neg_sig),1);
    tmp_text2=cellfun(@num2str,(num2cell(all_ps(neg_sig)')),'UniformOutput',0);
    sig_part{1}=strcat(strcat(tmp_text,tmp_text2)');
      end
    end
end

pos_check=0;
sig_part{2}={'no poscluster'};
pos_sig=[];

if isfield(stat, 'posclusters')
    if ~isempty(stat.posclusters)
      all_ps=[stat.posclusters.prob]
      pos_sig=find(all_ps<=alpha);
      if any(pos_sig)
    pos_check=1;
    tmp_text=repmat({'poscluster:p='},numel(pos_sig),1);
    tmp_text2=cellfun(@num2str,(num2cell(all_ps(pos_sig)')),'UniformOutput',0);
    sig_part{2}=strcat(strcat(tmp_text,tmp_text2)');
      end
    end
end

sig_text=strcat(sig_part{1},sig_part{2})