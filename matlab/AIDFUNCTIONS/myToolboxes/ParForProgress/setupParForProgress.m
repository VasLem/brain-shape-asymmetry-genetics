function [path,ID] = setupParForProgress(N)
         warning off; %#ok<*WNOFF>
         path = pwd;
         
         if strcmp(computer,'PCWIN64')
            t = 'C:/PARFORTMP/';
         elseif strcmp(computer,'GLNXA64')
            %cd /home/pclaes4/Documents/MATLAB/PARFORTMP/;
            t = (['/home/' getUserName '/Documents/MATLAB/PARFORTMP/']);
         else
            %cd /Users/pclaes4/Documents/MATLAB/PARFORTMP/;
            t = (['/Users/' getUserName '/Documents/MATLAB/PARFORTMP/']);
         end
         mkdir(t);
         cd(t);
         % generate random 
         ID = num2str(randi(10000,1));
         mkdir(ID);cd(ID);
         parfor_progress(N);
         warning on; %#ok<*WNON>
end