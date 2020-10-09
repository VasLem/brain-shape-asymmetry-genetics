function [f,TOTSNP] = METAManhattanplot(CHRRES,DB,PeakRS,index,cmap)
    f = figure;hold on;
    f.Color = [1 1 1];
    if isempty(cmap), cmap = colormap('lines'); end
    %cmap = colormap('colorcube');
  
%     cmap = cmap + repmat([0.5 0 0],size(cmap,1),1);
%     cmap(cmap>1) = 1;
%     
%     cmap = colormap('gray');
%     cmap = [cmap(1:3:end,:);cmap(end,:)];
%     cmap = cmap + repmat([0.5 0 0],size(cmap,1),1);
%     cmap(cmap>1) = 1;
    

%     step = 0.1;
%     new = cmap(1:7,:);
%     newmap = [new-0.3;new-0.1;new+0.1;new+0.3];
%     newmap(newmap>1) = 1;
%     newmap(newmap<0) = 0;
%     cmap = newmap;
% 
%     cmap = colormap('parula');
%     cmap = [cmap(1:3:end,:);cmap(end,:)];
    nrCHR = length(CHRRES);
    %if nargin<3, index = 1:63;end
    
    
%     index = 1:63;
%     PeakRS = [];
%     DB = 1;
    
    nMOD = length(index);
    X = 0;
    XTick = zeros(1,nrCHR);
    Xlab = cell(1,nrCHR);
    if ~isempty(PeakRS)
        peak = true;
    else
        peak = false;
    end
    maxp = 0;
    minp = -log10(1);
    TOTSNP = 0;
    for c=1:1:nrCHR
        % c=1;
        C = CHRRES{c};
        if peak, [ind12,ind21] = vlookupFast(PeakRS,C.RS);end
        C.nrSNP = length(C.POS);
        switch DB
            case -inf
                C.P = min(C.PFASTVAL,[],3)';
            case +inf
                C.P = max(C.PFASTVAL,[],3)';
            case 1000
                C.P = C.PFASTVAL';
            otherwise
                C.P = squeeze(C.PFASTVAL(:,:,DB))';
        end
        
        
        x = X+1:X+C.nrSNP;
        XTick(c) = mean(x);
        X = X+C.nrSNP;
        %indmaf = find(C.MAF>=0.01);
        indmaf = 1:C.nrSNP;
        TOTSNP = TOTSNP+C.nrSNP;
        if peak
           ind12 = intersect(ind12,indmaf);
        else
           ind12 = [];
        end
        for m=1:1:1
            %m=1
            y = -log10(C.P(:,index(m)));
            ind = find(y>minp); 
            ind = intersect(indmaf,ind);
            plot(x(ind),y(ind),'Color',cmap(c,:),'LineStyle','none','Marker','.');
            maxy = max(y);
            maxp = max(maxp,maxy);
        end
        if peak&&~isempty(ind12) 
           ps = min(C.P(ind12,:)');
           plot(x(ind12),-log10(ps),'ko','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5]);
        end

        plot(x,-log10(5*10^(-8))*ones(1,length(x)),'k:','LineWidth',1)
        %plot(x,-log10(7.7*10^(-8))*ones(1,length(x)),'k-','LineWidth',1)
        plot(x,-log10(2.8238*10^(-8))*ones(1,length(x)),'k-.','LineWidth',1)
        plot(x,-log10(1.2821*10^(-9))*ones(1,length(x)),'k-','LineWidth',1)
        %plot(x,-log10(5*10^(-8))*ones(1,length(x)),'k-.','LineWidth',1)
        if c<23
           Xlab{c} = num2str(c);
        else
           Xlab{c} = 'X';
        end
    end
    f.CurrentAxes.XTick = XTick;
    f.CurrentAxes.XTickLabel = Xlab;
    %f.CurrentAxes.YLim = [minp maxp+5];
    f.CurrentAxes.YLim = [minp 100];
    %f.CurrentAxes.YTick = [1 5 7 8 10:5:maxp+5];
    f.CurrentAxes.YTick = [5 8 10:5:25];
    %set(f.CurrentAxes,'YScale','log');
    f.Position = [1 41 1920 540];
    xlabel('Chromosome');
    ylabel('-log_1_0(p)');
end