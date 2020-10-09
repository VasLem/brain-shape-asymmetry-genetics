function [pc,lc,F] = decompComprimise(C)
    % 3.5 PCA of the compromise
    [pc,lc]=eigen(C);
    F=pc*diag(sqrt(lc));
    minF=min(F);
    maxF=max(F);
    nS = size(C,1);
    names = cell(1,nS);
    for i=1:1:nS
        names{i} = num2str(i);
    end
%     figure;clf
%     ax1=1;ax2=2;
%     plotxyrange(F,ax1,ax2,...
%        'Compromise as barycenter of Studies',...
%        char(names),...
%        [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)]);
%     figure;clf
%     ax1=1;ax2=3;
%     plotxyrange(F,ax1,ax2,...
%        'Compromise as barycenter of Studies',...
%        char(names),...
%        [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)]);
%     figure;clf
%     ax1=2;ax2=3;
%     plotxyrange(F,ax1,ax2,...
%        'Compromise as barycenter of Studies',...
%        char(names),...
%        [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)]);

end