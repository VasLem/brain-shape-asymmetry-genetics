function [path,ID] = setupParForProgress(N)
         warning off; %#ok<*WNOFF>
         path = pwd;
         if strcmp(computer,'PCWIN64')
            cd C:/PARFORTMP/;
         elseif strcmp(computer,'GLNXA64')
            %cd /home/pclaes4/Documents/MATLAB/PARFORTMP/;
            cd(['/home/' getUserName '/Documents/MATLAB/PARFORTMP/']);
         else
            %cd /Users/pclaes4/Documents/MATLAB/PARFORTMP/;
            cd(['/Users/' getUserName '/Documents/MATLAB/PARFORTMP/']);
         end
         % generate random 
         ID = num2str(randi(10000,1));
         mkdir(ID);cd(ID);
         parfor_progress(N);
         warning on; %#ok<*WNON>
end