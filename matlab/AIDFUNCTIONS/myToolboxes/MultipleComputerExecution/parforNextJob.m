function todo = parforNextJob(jobs,path)
         if nargin==2, cd(path); end
         files = dir('*.mat');    
         if ~isempty(files)
            nrfiles = length(files);
            done = zeros(1,nrfiles);
            for i=1:1:nrfiles
                done(i) = str2double(files(i).name(1:end-6));
            end
         else
            done = [];
         end
         todo = setdiff(jobs,done);
         
         
         files = {files(~[files.isdir]).name};







end