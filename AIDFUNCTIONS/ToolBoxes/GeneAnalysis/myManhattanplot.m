function f = myManhattanplot(TEST,CHRID)
    f = figure;hold on;
    str = {'r.' 'k.' 'b.' 'g.' 'm.'};
    ax = {};
    indT = [1 3 5 7 2 4 6];
    models = {'aa aAAA' 'aaaA AA' 'aaAA aA' 'aa aA AA' 'aa AA' 'aa aA' 'aA AA'}; 
    for t=1:1:7
        ax{t} = subplot(4,2,indT(t));hold on;grid on;title(models{t});view(ax{t},0,0);
    end
    counter = zeros(1,7);
    chr = unique(CHRID);
    nrCHR = length(chr);
    LC = size(TEST,1);
    for c=1:nrCHR
        %c=1;
        %index = find(CHRID==chr(c));
        chrindex = find(CHRID==chr(c));
        tmp = TEST(:,:,chrindex);
        p = -log10(tmp);
        strc = str{mod(c,5)+1};
        nSNP = size(p,3);
        for t=1:1:7
            % t=1;
            tp = squeeze(p(:,t,:));
            x = counter(t)+1:counter(t)+nSNP;
            tax = ax{t};
            for i=1:1:LC;
                %i=1;
                y = (i-1)*10*ones(1,nSNP);
                z = tp(i,:);
                index = find(z>-log10(0.01));
                plot3(tax,x(index),y(index),z(index),strc,'MarkerSize',10);
            end
            counter(t) = counter(t)+nSNP;
            x = 1:counter(t);
            %plot3(tax,x,zeros(1,length(x)),-log10(10^(-6))*ones(1,length(x)),'c--','LineWidth',1);
            plot3(tax,x,zeros(1,length(x)),-log10(5*10^(-8))*ones(1,length(x)),'c-.','LineWidth',1);
            drawnow;
        end
            %pTOT = single(cat(3,pTOT,p));
    end
    for t=1:1:7
        set(ax{t},'zlim',[3 15]);
    end
end