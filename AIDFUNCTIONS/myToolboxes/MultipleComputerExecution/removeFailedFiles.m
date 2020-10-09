function removeFailedFiles(path)
    if nargin<1, path = pwd; end
    cd(path);
    COMPID = getComputerID;
    files = dir(['*_' num2str(COMPID) '_2.mat']);
    for i=1:1:length(files)
        if files(i).bytes==0, delete(files(i).name); disp([files(i).name ' Deleted']); end
    end
end