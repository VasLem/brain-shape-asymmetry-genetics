function [f] = BRAINManhattanplot(POS,CHR,P,pGW,pSW,mode,cmap,peaks)
    f = figure;hold on;f.Color = [1 1 1];
    if nargin<8, peaks = []; end
    if nargin<7, cmap = colormap('lines'); end
    if nargin<6, mode = 'best'; end
    if nargin<5,pSW = 1e-10; end
    if nargin<4, pGW = 5e-8;end
    
    if pSW<1, pSW = -log10(pSW);end
    if pGW<1, pGW = -log10(pGW);end
    
    MAX = max(P(:));
    MAX = min(MAX,80);
    f.CurrentAxes.YLim = [0 MAX];
    f.CurrentAxes.YTick = [0 3 5 8 12 20 30 40 50 60 70];
    f.CurrentAxes.YTickLabel = abs(f.CurrentAxes.YTick);
    
    nSNP = length(POS);X = 1:nSNP;
    CHRID = unique(CHR);
    nCHR = length(CHRID);
    XTick = zeros(1,nCHR);
    Xlab = cell(1,nCHR);
    for c=1:1:nCHR
       %c=1;
       index = find(CHR==CHRID(c)); 
       XTick(c) = median(X(index)); 
       if c<23
          Xlab{c} = num2str(c);
       else
          Xlab{c} = 'X';
       end
    end
    plot(X,ones(1,nSNP)*pGW,'k:','LineWidth',1);
    plot(X,ones(1,nSNP)*pSW,'k-','LineWidth',1);
    f.CurrentAxes.XLim = [1 nSNP];
    f.CurrentAxes.XTick = XTick;
    f.CurrentAxes.XTickLabel = Xlab;
    f.Position = [1 41 2500 1000];
    xlabel('Chromosome');
    ylabel('-log_10(p)');
    colormap(f,cmap);
    switch lower(mode)
        case 'best'
            Pmax = max(P,[],2);
            scatter(X,Pmax,5,CHR,'filled');
        otherwise
            Pmax = max(P,[],2);
            scatter(X,Pmax,5,CHR,'filled');
            for m=1:size(P,2)
               pvalues = P(:,m);
               index = find(pvalues>3);
               scatter(X(index),pvalues(index),5,CHR(index),'filled'); 
            end
    end
    if isempty(peaks), return; end
    [ind12,ind21] = vlookupFast([uint32(peaks.CHR) uint32(peaks.POS)],[uint32(CHR) uint32(POS)],true);
    x = X(ind12);
    y = P(ind12,:);
    y = max(y,[],2);
    y(y>80) = 80;
    plot(x,y,'k.','MarkerSize',20);
    
end