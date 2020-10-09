function parforLocking(lockpath)
         cd(lockpath);
         lock = true;
         hold = true;
         while hold
            files = dir('LOCKED*'); % dealing with conflicted copies
            if isempty(files)
               save('LOCKED','lock');
               hold = false;
            else
               pause(10*rand(1));
            end
         end
end