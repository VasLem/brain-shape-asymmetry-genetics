function [path,ID] = setupParForProgress(N)
         warning off; %#ok<*WNOFF>
         path = pwd;
         if ~startsWith(getComputerName, 'mic')
             if strcmp(computer,'PCWIN64')
                direc = 'C:/PARFORTMP/';
             elseif strcmp(computer,'GLNXA64')
                %cd /home/pclaes4/Documents/MATLAB/PARFORTMP/;
                direc = ['/home/' getUserName '/Documents/MATLAB/PARFORTMP/'];
             else
                %cd /Users/pclaes4/Documents/MATLAB/PARFORTMP/;
                direc = ['/Users/' getUserName '/Documents/MATLAB/PARFORTMP/'];
             end
         else
             direc = ['/usr/local/avalok/tmp/' getUserName '/MATLAB/PARFORTMP/'];
         end
         mkdir(direc);
         cd(direc);
         % generate random
         ID = num2str(randi(10000,1));
         mkdir(ID);cd(ID);
         parfor_progress(N);
         warning on; %#ok<*WNON>
end