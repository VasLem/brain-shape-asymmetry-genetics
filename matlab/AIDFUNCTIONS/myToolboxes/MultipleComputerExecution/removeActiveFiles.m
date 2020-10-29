function removeActiveFiles(path)
    if nargin<1, path = [pwd '/']; end
    delete([path '*_' num2str(getComputerID) '_0.mat']);
end