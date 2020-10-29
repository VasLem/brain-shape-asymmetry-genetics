function [f] = METAMiamiplotPerCHR(C,PeakInfo,index,color,pGW,pSW,ID)
    f = figure;hold on;
    f.Color = [1 1 1];
    if isempty(color), color = colormap('lines'); end
    nMOD = length(index);
    nSNP = length(C.POS);
    
    falsecolor = [62 38 168]/255;
    intermediatecolor = [248 186 61]/255;
    truecolor = [249 251 20]/255;
    
    x = 1:nSNP;
    %x = C.POS(:,2)';
    f.CurrentAxes.XLim = [x(1) x(end)];
    %XTick = x(1):(x(end)-x(1))/10000:x(end);
    XTick = 1:10000:nSNP;
    Xlab = cell(1,length(XTick));
    for i=1:1:length(XTick)
        Xlab{i} = num2str(C.POS(XTick(i),2));
        %Xlab{i} = num2str(XTick(i));
    end
    if ~isempty(PeakInfo)
        peak = true;
    else
        peak = false;
    end
    if peak 
        %index = find(PeakInfo.CHR==ID);
        [ind12,ind21] = vlookupFast(PeakInfo.RS,C.RS(:,2));
    end
    P1 = squeeze(C.PFASTVAL(:,:,1))';
    P2 = squeeze(C.PFASTVAL(:,:,2))';
    maxp = 0;
    minp = -log10(1);
    for m=1:1:nMOD
        %m=1
        y = -log10(P1(:,index(m)));
        ind = find(y>minp); 
        plot(x(ind),y(ind),'Color',color,'LineStyle','none','Marker','.');
        maxy = max(abs(y));
        maxp = max(maxp,maxy);
        y = log10(P2(:,index(m)));
        ind = find(abs(y)>minp); 
        plot(x(ind),y(ind),'Color',color,'LineStyle','none','Marker','.');
        maxy = max(abs(y));
        maxp = max(maxp,maxy);
    end
    if peak&&~isempty(ind12)
        for p=1:1:length(ind12)
            % p=1;
            xp = x(ind12(p));
            switch PeakInfo.BESTSym(ind21(p))
                case 0
                    switch PeakInfo.SUGSym(ind21(p))>0 
                        case true
                            pcolor = intermediatecolor;
                            showrs = true;
                        case false
                            pcolor = falsecolor;
                            showrs = false;
                    end
                case 1
                    pcolor = truecolor;
                    showrs = true;
            end
            switch PeakInfo.BestDB(ind21(p))
                case 1
                    yp = -log10(PeakInfo.topP(ind21(p)));
                case 2
                    yp = log10(PeakInfo.topP(ind21(p)));
            end
            plot(xp,yp,'ko','MarkerSize',5,'MarkerFaceColor',pcolor);
            switch mod(p,2)
                case 0
                    HA = 'right';
                case 1
                    HA = 'left';
            end
            if showrs, text(xp,yp,PeakInfo.RS{ind21(p)},'FontSize',8,'HorizontalAlignment',HA);end
        end
    end
    f.CurrentAxes.YLim = [-1*maxp maxp];
    X = 1:1000:nSNP;
    plot(X,-log10(pGW)*ones(1,length(X)),'k:','LineWidth',1)
    plot(X,log10(pGW)*ones(1,length(X)),'k:','LineWidth',1)
    plot(X,-log10(pSW)*ones(1,length(X)),'k--','LineWidth',1)
    plot(X,log10(pSW)*ones(1,length(X)),'k--','LineWidth',1)
    plot(X,0*ones(1,length(X)),'k-','LineWidth',1)    
    
    
    f.CurrentAxes.XTick = XTick;
    f.CurrentAxes.XTickLabel = Xlab;
    f.CurrentAxes.XTickLabelRotation = 90;
    f.CurrentAxes.YTickLabel = abs(f.CurrentAxes.YTick);
    
    %set(f.CurrentAxes,'YScale','log');
    f.Position = [508         203        2595         869];
    xlabel('Base-Pair position');
    ylabel('-log_1_0(p)');
    title(['CHROMOSOME ' num2str(ID)]);
end