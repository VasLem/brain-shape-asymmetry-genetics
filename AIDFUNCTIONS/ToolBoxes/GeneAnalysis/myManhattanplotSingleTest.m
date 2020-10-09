function f = myManhattanplotSingleTest(TEST,CHRID,pCrit)
    if nargin<3, pCrit = 5*10^-8;end
    f = figure;hold on;view(0,0);
    str = {'r.' 'k.' 'b.' 'g.' 'm.'};
    models = {'aa aA AA'}; 
    title(models{1});
    ax = gca;
    chr = unique(CHRID);
    nrCHR = length(chr);
    LC = size(TEST,1);
    counter = 0;
    for c=1:nrCHR
        %c=1;
        %index = find(CHRID==chr(c));
        chrindex = find(CHRID==chr(c));
        tmp = TEST(:,chrindex);
        tp = -log10(tmp);
        strc = str{mod(c,5)+1};
        nSNP = size(tp,2);
        x = counter+1:counter+nSNP;
        for i=1:1:LC;
                %i=1;
                y = (i-1)*10*ones(1,nSNP);
                z = tp(i,:);
                index = find(z>-log10(1));
                plot3(ax,x(index),y(index),z(index),strc,'MarkerSize',10);
        end
        counter = counter+nSNP;
        x = 1:counter;
        %plot3(tax,x,zeros(1,length(x)),-log10(10^(-6))*ones(1,length(x)),'c--','LineWidth',1);
        plot3(ax,x,zeros(1,length(x)),-log10(pCrit)*ones(1,length(x)),'c-.','LineWidth',1);
        drawnow;
    end
    %set(ax,'zlim',[1 20]);
end