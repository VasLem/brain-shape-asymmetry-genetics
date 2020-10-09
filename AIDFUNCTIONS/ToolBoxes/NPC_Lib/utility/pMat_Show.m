function []=pMat_Show(p,alpha,labels,nfig,rows_cols)
%%%%%%%%%%%%%%%%
%pMat_Show(p,alpha,labels,nfig,rows_cols)
%produce the coroled map from the matrix p.
%
%INPUTS
%p= (B+1)xMYxMXxNS matrix of p-values.
%   (B+1): number of random permuations (the function just make use of the
%   last one, therefore you can simply iNPut a 1xMYxMXxNS matrix.
%   MY   : number of Y variables
%   MX   : number of X variables
%   MS   : number of Stratas
%Use signed p-values for test with different tails.
%alpha    = significance level (alpha=.05 by default).
%labels   = see label_std (not required iNPut argument).
%nfig     = number of figure to open.
%rows_cols= diplays the subplots in rows_cols(1) rows and rows_cols(2)
%           columns (not required iNPut argument)
% Livio Finos (2006)


if nargin<3
    labels=label_std(p(end,:,:,:));    
else
    if isempty(labels)
        labels=label_std(p(end,:,:,:));   
    else
        if not(isfield(labels, 'dims')) ...
                | not(isfield(labels, 'dimslabel'))
            if isfield(labels, 'labels')
            labels=labels.labels;
            else
                labels=label_std(p(end,:,:,:));   
            end
        end
    end
end

if not(isfield(labels, 'title'))
    if length(labels.dims)>3
        labels.title=labels.dims{4};
    else
    for i=1:size(p,3)
        labels.title(i)={['DIM' num2str(i)] };
    end
    end
end

 p(abs(p)>alpha)=sign(p(abs(p)>alpha)).*alpha*(1+6/70);
% p(abs(p)<10^-16)=sign(p(abs(p)<10^-16))*10^-16;
 
 if nargin<5
    sq=ceil(sqrt(size(p,4)));
    rows=ceil(size(p,4)./sq);
 else
     sq=rows_cols(1);
     rows=rows_cols(2);
 end

  
 if nargin>=4
     figure(nfig) %'PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
 else
     figure
 end
 
 
load pMat_colbar
if (sum(p(end,:)>0)*sum(p(end,:)<0))>0
    MAP=P_bi;
    range=[-alpha*(1+7/70) alpha*(1+7/70)];
    Ticks=[-alpha*(1+ 6/70) -alpha:alpha/4:alpha alpha*(1+ 6/70)];
    TicksLabel={'NS' num2str(alpha) num2str(alpha*3/4) num2str(alpha*2/4) num2str(alpha/4) num2str(0) num2str(alpha/4) num2str(alpha*2/4) num2str(alpha*3/4) num2str(alpha) 'NS' };
else
    MAP=P_hot;
    range=[0 alpha*(1+7/70)];
    Ticks=[0 alpha/4:alpha/4:alpha alpha*(1+ 6/70)];
    TicksLabel={num2str(0) num2str(alpha/4) num2str(alpha*2/4) num2str(alpha*3/4) num2str(alpha) 'NS'};
end

     
for i=1:size(p,4)
    temp=squeeze(p(end,:,:,i));
    if size(temp,1)==1
        temp=temp';
    end
    subplot(sq,rows,i),imagesc(temp,range);  colormap(MAP);
    if size(p,4)>1,title([labels.title{i}]);end
      if mod(i,rows)==1, 
          ylabel(labels.dimslabel{2});
      end
      if i>size(p,3)-rows
          xlabel(labels.dimslabel{3});
      end
   set(gca,'Box','off')
   set(gca,'TickDir','out')
    set(gca,'FontSize',min(ceil(240/length(labels.dims{2})/1.5),15))
    set(gca,'Xtick',1:size(p,3))
    set(gca,'Ytick',1:size(p,2))
    set(gca,'XTickLabel',labels.dims{3}(get(gca,'Xtick')))
%    set(gca,'Rotation',90)

    set(gca,'YTickLabel',labels.dims{2}(get(gca,'Ytick')))
    
% for iii=1:length(labels.dims{2})
%     text(iii,length(labels.dims{3})+1,labels.dims{2}(iii),'rotation',90,'HorizontalAlignment','right')
% end
% 
    if mod(i,rows)==0, 
    
    colorbar('YTick',Ticks,'YTickLabel', TicksLabel);
    
    end
    grid on
end 
tilefigs
