function plotExp(resultList, parameters, reason)
pGroups = ["permFF","permIF","permDF"; "FF","IF", "DF"; "F", "I", "D"]'; 
for p=1:3
    for i=1:size(resultList,2)
        for j=1:3
            toPlot(i,j) = resultList(i).Total.(pGroups(j, p));
        end
    end
    figure();
    plot(parameters, toPlot(:, 1),'r');
    hold on;
    plot(parameters, toPlot(:, 2), 'g');
    plot(parameters, toPlot(:, 3),'b');
    xlabel(reason);
    legend(pGroups(:,p));
    title(["ANOVA2-way asymmetry significance~"  reason]);
    hold off;
end

end