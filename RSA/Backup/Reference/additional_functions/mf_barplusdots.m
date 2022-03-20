% plotting function for plotting bars plus data as single points
% data: organized in sub x conds

function mf_barplusdots(data,barlabels,ylabel_str)
 bar(mean(data))
 hold on
 % for scatter plot
 % xdata: 1:number of bars, with jitter for overlapping data
 xdata=[];
for  col=1:size(data,2)
xtmp=ones(size(data,1),1).*col;
ytmp=data(:,col);
[C,ia,ic] = unique(ytmp);
non_unique=setdiff(ic,ia);
    if ~isempty(non_unique)
        for i=1:numel(non_unique)
            same_ind=ytmp==ytmp(non_unique(i));
            num_sam=sum(same_ind);
            even=(floor(num_sam/2)==(num_sam/2));
            jit=(floor(num_sam/2).*(0.025))-(even.*(0.025/2));
            newx=-jit:0.025:jit;
            xtmp(same_ind)=xtmp(same_ind)+newx';
            clear newx
        end
    end
    xdata=[xdata;xtmp];
end    
 %[ones(18,1);ones(18,1).*2;ones(18,1).*3;ones(18,1).*4] 
 
 % ydata: data values
 ydata=reshape(data,[],1);
 scatter(xdata,ydata,5,'k','filled');
 xticklabels(barlabels)
 ylabel(ylabel_str)