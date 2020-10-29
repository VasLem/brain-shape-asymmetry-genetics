function removeActiveFilesv2(varargin)
    Input = find(strcmpi(varargin, 'path'));
    if isempty(Input),path = [pwd '/'];else path = varargin{Input+1};end
    Input = find(strcmpi(varargin, 'id'));
    if isempty(Input),id = getComputerID;else id = varargin{Input+1};end
    if ~strcmp(path(end),'/'),path = [path '/'];end
    delete([path '*_' num2str(id) '_0.mat']);
end