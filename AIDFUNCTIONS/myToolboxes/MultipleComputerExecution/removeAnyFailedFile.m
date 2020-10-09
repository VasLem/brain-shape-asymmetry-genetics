function removeAnyFailedFile(path)
    if nargin<1, path = pwd; end
    cd(path);
    files = dir('*_2.mat');
    for i=1:1:length(files)
        delete(files(i).name); disp([files(i).name ' Deleted']);
    end
end