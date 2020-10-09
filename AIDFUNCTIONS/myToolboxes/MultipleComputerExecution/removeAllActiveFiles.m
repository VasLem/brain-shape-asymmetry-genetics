function removeAllActiveFiles(path)
    if nargin<1, path = [pwd '/']; end
%     bkpath = pwd;
%     cd(path);
%     files = dir(['*_' num2str(COMPID) '_0.mat']);
%     for i=1:1:length(files)
%         if files(i).bytes==0, delete(files(i).name); disp([files(i).name ' Deleted']); end
%     end
%     cd(bkpath);
    delete([path '*_0.mat']);
end