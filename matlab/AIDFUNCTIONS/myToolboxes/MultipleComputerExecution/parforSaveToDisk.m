function parforSaveToDisk(path,filename,in)
         if ~strcmp(path(end),'/')
            save([path '/' filename],'in','-v7.3');
         else
            save([path filename],'in','-v7.3'); 
         end
end