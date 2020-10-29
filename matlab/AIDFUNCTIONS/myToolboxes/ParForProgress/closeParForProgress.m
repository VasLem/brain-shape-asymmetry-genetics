function closeParForProgress(path,ID)
         try
             if strcmp(computer,'PCWIN64')
                cd C:/PARFORTMP/;
             elseif strcmp(computer,'GLNXA64')
                %cd /home/pclaes4/Documents/MATLAB/PARFORTMP/;
                cd(['/home/' getUserName '/Documents/MATLAB/PARFORTMP/']);
             else
                %cd /Users/pclaes4/Documents/MATLAB/PARFORTMP/;
                cd(['/Users/' getUserName '/Documents/MATLAB/PARFORTMP/']);
             end
         catch
             cd(['/usr/local/avalok/tmp/' getUserName '/MATLAB/PARFORTMP/']);
         end
         cd(ID);
         parfor_progress(0);
         %cd ..
         %rmdir(ID);
         cd(path);
end