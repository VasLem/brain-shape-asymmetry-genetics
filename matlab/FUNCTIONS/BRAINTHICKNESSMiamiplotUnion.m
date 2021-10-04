function [f] = BRAINTHICKNESSMiamiplotUnion(BPOS,BCHR,BP,TPOS,TCHR,TP,pGW,pSWB,pSWT,mode,cmap)
    
    if pSWB<1, pSWB = -log10(pSWB);end
    if pSWT<1, pSWT = -log10(pSWT);end
    if pGW<1, pGW = -log10(pGW);end
    
    
    pSW = cell(1,2);
    pSW{1} = pSWB;
    pSW{2} = pSWT;
    pSUG = -log10(5e-7);
    
    CHRID = unique(TCHR);
    nCHR = length(CHRID);
    CHR = [];
    POS = [];
    [path,ID] = setupParForProgress(nCHR);tic;
    parfor c=1:nCHR
        % c=1;
       index=find(TCHR==CHRID(c));
       tmpposT = TPOS(index);
       index=find(BCHR==CHRID(c));
       tmpposB = BPOS(index);
       tmppos = union(tmpposT,tmpposB); 
      
       POS = [POS;uint32(tmppos(:))];
       CHR = [CHR;uint8(double(CHRID(c))*ones(length(tmppos),1))];
      parfor_progress;  
    end
    closeParForProgress(path,ID);toc;
   
    PP = cell(1,2);
    PP{2} = zeros(size(CHR,1),size(TP,2));
    [ind12,ind21] = vlookupFast([uint32(TCHR) TPOS],[uint32(CHR) POS],true);
    PP{2}(ind12,:) = TP(ind21,:);
   
    PP{1} = zeros(size(CHR,1),size(BP,2));
    [ind12,ind21] = vlookupFast([uint32(BCHR) BPOS],[uint32(CHR) POS],true);
    PP{1}(ind12,:) = BP(ind21,:);
   
    
    % creating figure
    f = figure;hold on;f.Color = [1 1 1];f.Position = [1 41 2500 1000];
    MAX = [70 70];
    f.CurrentAxes.YLim = [-1*MAX(1) MAX(2)];
    f.CurrentAxes.YTick = [-70 -60 -50 -40 -30 -20 -10 -4 4 10 20 30 40 50 60 70];
    f.CurrentAxes.YTickLabel = abs(f.CurrentAxes.YTick);
    
    nSNP = length(POS);X = 1:nSNP;
    
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
                for m=1:size(P,2)
                   pvalues = P(:,m);
                   index = find(pvalues>=5);
                   scatter(X(index),signs*pvalues(index),5,CHR(index),'filled'); 
                end
        end
    end
    
    plot(X,ones(1,nSNP)*pGW,'k--','LineWidth',1);
    %plot(X,ones(1,nSNP)*pSUG,'k:','LineWidth',1);
    plot(X,ones(1,nSNP)*pSW{1},'k-','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pGW,'k--','LineWidth',1);
    %plot(X,-1*ones(1,nSNP)*pSUG,'k:','LineWidth',1);
    plot(X,-1*ones(1,nSNP)*pSW{2},'k-','LineWidth',1);
    plot(X,zeros(1,nSNP),'k-','LineWidth',1);
    
end