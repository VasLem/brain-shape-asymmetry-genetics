function out = dcRetrieveInputFolders(path)
         outstr = genpath(path);
         index = strfind(outstr,':');
         n = length(index);
         out = cell(1,length(index));
         for i=1:1:n
             if i==1
                out{i} = outstr(1:index(i)-1);
             else
                out{i} = outstr(index(i-1)+1:index(i)-1);
             end
             if ~strcmp(out{i}(end),'/'), out{i} = [out{i} '/'];end
         end
end