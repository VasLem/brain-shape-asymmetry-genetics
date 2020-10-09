function [f] = BRAINFACEMiamiplot(BPOS,BCHR,BP,FPOS,FCHR,FP,pGW,pSWB,pSWF,mode,cmap,peaks,peaks2)
    
    if pSWB<1, pSWB = -log10(pSWB);end
    if pSWF<1, pSWF = -log10(pSWF);end
    if pGW<1, pGW = -log10(pGW);end
    
    
    pSW = cell(1,2);
    pSW{1} = pSWB;
    pSW{2} = pSWF;
    pSUG = -log10(5e-7);
    
    
    % lookfor interesection
    [ind12,ind21] = vlookupFast([uint32(FCHR) FPOS],[uint32(BCHR) BPOS],true);
    CHR = FCHR(ind21);
    POS = FPOS(ind21);
    PP = cell(1,2);
    PP{1} = BP(ind12,:);
    PP{2} = FP(ind21,:);
%     PP{1} = single(BP(ind12,:))/10000;
%     PP{2} = single(FP(ind21,:))/10000;
    
    % creating figure
    f = figure;hold on;f.Color = [1 1 1];f.Position = [1 41 2500 1000];
    MAX = [70 70];
    f.CurrentAxes.YLim = [-1*MAX(1) MAX(2)];
    f.CurrentAxes.YTick = [-70 -60 -50 -40 -30 -20 -10 -5 5 10 20 30 40 50 60 70];
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
    plot(X,ones(1,nSNP)*pSUG,'k:','LineWidth',1);
    plot(X,ones(1,nSNP)*pSW{1},'k-','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pGW,'k:','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pSUG,'k:','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pSW{2},'k-','LineWidth',1);
    plot(X,zeros(1,nSNP),'k-','LineWidth',1);
    f.CurrentAxes.XLim = [1 nSNP];
    f.CurrentAxes.XTick = XTick;
    f.CurrentAxes.XTickLabel = Xlab;

    xlabel('Chromosome');
    ylabel('-log_1_0(p)');
    colormap(f,cmap);
    for i=1:2
        % i=1;
        switch i
            case 1
                P = PP{1};
                signs = 1;
            case 2
                P = PP{2};
                signs = -1;
        end
        switch lower(mode)
            case 'best'
%                 Pmax = max(P,[],2);
%                 scatter(X,signs*Pmax,5,CHR,'filled');
            otherwise
                %Pmax = max(P,[],2);
                %scatter(X,signs*Pmax,5,CHR,'filled');
%                 for m=1:size(P,2)
%                    pvalues = P(:,m);
%                    index = find(pvalues>=3);
%                    scatter(X(index),signs*pvalues(index),5,CHR(index),'filled'); 
%                 end
        end
    end
    if isempty(peaks), return; end
    for i=1:length(peaks.RS)
        % i=1;
        chr = peaks.CHR(i);
        pos = double(peaks.POS(i));
        index = find(CHR==chr);
        dist = abs(double(POS(index))-pos);
        [~,minind] = min(dist);
        x = X(index(minind));
        y = 0;
        plot(x,y,'k.','MarkerSize',15);
    end
    if isempty(peaks2), return; end
    for i=1:length(peaks2.RS)
        % i=1;
        chr = peaks2.CHR(i);
        pos = double(peaks2.POS(i));
        index = find(CHR==chr);
        dist = abs(double(POS(index))-pos);
        [~,minind] = min(dist);
        x = X(index(minind));
        y = 0;
        plot(x,y,'k.','MarkerSize',15);
    end
end