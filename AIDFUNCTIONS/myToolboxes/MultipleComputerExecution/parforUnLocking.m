function parforUnLocking(lockpath)
         cd(lockpath);
         files = dir('LOCKED*');% dealing with conflicted copies through dropbox
         if isempty(files), return; end
         for i=1:1:length(files)
             try
                 delete(files(i).name);
             catch
             end
         end
end