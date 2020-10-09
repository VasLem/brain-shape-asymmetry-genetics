function [f] = BRAINMiamiplot(POS,CHR,PLH,PRH,pGW,pSW,mode,cmap)
    f = figure;hold on;f.Color = [1 1 1];
    if nargin<8, cmap = colormap('lines'); end
    if nargin<7, mode = 'best'; end
    if nargin<6,pSW = 1e-10; end
    if nargin<5, pGW = 5e-8;end
    
    if pSW<1, pSW = -log10(pSW);end
    if pGW<1, pGW = -log10(pGW);end
    
    MAX = max([max(PLH(:)) max(PRH(:))]);
    MAX = min(MAX,80);
    f.CurrentAxes.YLim = [-1*MAX MAX];
    f.CurrentAxes.YTick = [-70 -60 -50 -40 -30 -20 -12 -8 8 12 20 30 40 50 60 70];
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
    plot(X,-1*ones(1,nSNP)*pGW,'k:','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pSW,'k-','LineWidth',1);
    plot(X,zeros(1,nSNP),'k-','LineWidth',1);
    f.CurrentAxes.XLim = [1 nSNP];
    f.CurrentAxes.XTick = XTick;
    f.CurrentAxes.XTickLabel = Xlab;
    f.Position = [1 41 2500 1000];
    xlabel('Chromosome');
    ylabel('-log_1_0(p)');
    colormap(f,cmap);
    for i=1:2
        % i=1;
        switch i
            case 1
                P = PLH;
                signs = 1;
            case 2
                P = PRH;
                signs = -1;
        end
        switch lower(mode)
            case 'best'
                Pmax = max(P,[],2);
                scatter(X,signs*Pmax,5,CHR,'filled');
            otherwise
                Pmax = max(P,[],2);
                scatter(X,signs*Pmax,5,CHR,'filled');
                for m=1:size(P,2)
                   pvalues = P(:,m);
                   index = find(pvalues>6);
                   scatter(X(index),signs*pvalues(index),5,CHR(index),'filled'); 
                end
        end
    end
end