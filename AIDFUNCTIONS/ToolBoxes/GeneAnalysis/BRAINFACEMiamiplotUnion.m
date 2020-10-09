function [f] = BRAINFACEMiamiplotUnion(BPOS,BCHR,BP,FPOS,FCHR,FP,pGW,pSWB,pSWF,mode,cmap,peaks)
    
    if pSWB<1, pSWB = -log10(pSWB);end
    if pSWF<1, pSWF = -log10(pSWF);end
    if pGW<1, pGW = -log10(pGW);end
    
    
    pSW = cell(1,2);
    pSW{1} = pSWB;
    pSW{2} = pSWF;
    pSUG = -log10(5e-7);
    
    CHRID = unique(FCHR);
    nCHR = length(CHRID);
    CHR = [];
    POS = [];
    for c=1:nCHR
        % c=1;
       index=find(FCHR==CHRID(c));
       tmpposF = FPOS(index);
       index=find(BCHR==CHRID(c));
       tmpposB = BPOS(index);
       tmppos = union(tmpposF,tmpposB); 
      
       POS = [POS;uint32(tmppos(:))];
       CHR = [CHR;uint8(double(CHRID(c))*ones(length(tmppos),1))];
        
    end
    
    
%     
%     
%     % create union
%     un = union([uint32(FCHR) FPOS],[uint32(BCHR) BPOS],'rows');
%     CHR = uint8(un(:,1));
%     POS = un(:,1);
    
    PP = cell(1,2);
    PP{2} = zeros(size(CHR,1),size(FP,2));
    [ind12,ind21] = vlookupFast([uint32(FCHR) FPOS],[uint32(CHR) POS],true);
    PP{2}(ind12,:) = FP(ind21,:);
   
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
    
    if isempty(peaks), return; end
    up = 1;left = 1;
    for i=1:length(peaks.RS)
        %if peaks.TEST(i)==2, continue; end
        % i=3;
        chr = peaks.CHR(i);
        pos = double(peaks.POS(i));
        index = find(CHR==chr);
        dist = abs(double(POS(index))-pos);
        [~,minind] = min(dist);
        x = X(index(minind));
        MS = 8;
        markerline = ':';
        y = peaks.Brain.Proxy(i);
        if y>6.3
            plot([x;x],[0;y],markerline,'Color',cmap(chr,:),'LineWidth',1.5);
        else
            plot([x;x],[0;y],markerline,'Color',cmap(chr,:),'LineWidth',1.5);
        end
        y = peaks.Face.Proxy(i);
        if y>7.3
            plot([x;x],[0;-y],markerline,'Color',cmap(chr,:),'LineWidth',1.5);
        else
            plot([x;x],[0;-y],markerline,'Color',cmap(chr,:),'LineWidth',1.5);
        end
        y = peaks.Brain.Proxy(i);
        y = min([y 70]);
        if y>7.3
            plot(x,y,'o','Color',[0 0 0],'MarkerSize',MS,'LineWidth',2);
        else
            plot(x,y,'o','Color',[0 0 0],'MarkerSize',MS,'LineWidth',2);
        end
        y = peaks.Face.Proxy(i);
        y = min([y 70]);
        if y>7.3
           plot(x,-y,'o','Color',[0 0 0],'MarkerSize',MS,'LineWidth',2);
        else
           plot(x,-y,'o','Color',[0 0 0],'MarkerSize',MS,'LineWidth',2);
        end  
%         switch peaks.TEST(i)
%             case 2
%                 plot(x,0,'o','Color',[0 0 0],'MarkerSize',MS);
%             case 1
%                 plot(x,0,'o','Color',[0 0 0],'MarkerSize',MS);
%             otherwise
%                 plot(x,0,'o','Color',[0 0 0],'MarkerSize',MS);
%         end
        if isempty(peaks.LABELS{i}), continue; end
        y = up*3;up = up*-1;
        %plot([x;x],[0;.5*y],'k-','LineWidth',1);
        switch left
            case 1
                text(x,.55*y,peaks.LABELS{i},'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
            case -1
                text(x,.55*y,peaks.LABELS{i},'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');   
        end
        %plot([x;x],[0;y],'k-','LineWidth',1);
    end
%     for i=1:length(peaks.RS)
%         if peaks.TEST(i)==1, continue; end
%         % i=1;
%         chr = peaks.CHR(i);
%         pos = double(peaks.POS(i));
%         index = find(CHR==chr);
%         dist = abs(double(POS(index))-pos);
%         [~,minind] = min(dist);
%         x = X(index(minind));
%         y = 0;
%         switch peaks.TEST(i)
%             case 2
%                 plot(x,y,'.','Color',[0 0 0],'MarkerSize',20);
%             case 1
%                 plot(x,y,'.','Color',[0.5 0.5 0.5],'MarkerSize',20);
%             otherwise
%                 plot(x,y,'.','Color',[0.8 0.8 0.8],'MarkerSize',20);
%         end
%     end
    

    
    
%     if isempty(peaks2), return; end
%     for i=1:length(peaks2.RS)
%         % i=1;
%         chr = peaks2.CHR(i);
%         pos = double(peaks2.POS(i));
%         index = find(CHR==chr);
%         dist = abs(double(POS(index))-pos);
%         [~,minind] = min(dist);
%         x = X(index(minind));
%         y = 0;
%         plot(x,y,'k.','MarkerSize',15);
%     end
end