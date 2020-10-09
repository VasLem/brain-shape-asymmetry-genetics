function [DATA,time] = mapDataGroup(DATA,RMapper,NRMapper,followprogress)
    if nargin<4, followprogress = false; end
    N = length(DATA);
    if followprogress, [path,ID] = setupParForProgress(N);end
    tic;
    parfor i=1:N
        if ~isempty(RMapper), forRMapper = clone(RMapper); else, forRMapper = []; end
        if ~isempty(NRMapper), forNRMapper = clone(NRMapper); else, forNRMapper = []; end
        out = mapFace(DATA{i}.Shape,DATA{i}.InitTemplate,forRMapper,forNRMapper,false);
        DATA{i}.MappedShape = clone(out);
        if followprogress, parfor_progress;end
    end
    time = toc;
    if followprogress, closeParForProgress(path,ID);end
%     disp(['Group mapping done in ' num2str(time) ' seconds.']);
end