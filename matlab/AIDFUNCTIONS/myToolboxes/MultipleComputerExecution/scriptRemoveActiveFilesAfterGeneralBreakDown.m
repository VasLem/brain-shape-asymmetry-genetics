%% A script to remove all active files of computers in my list
path = '/home/pclaes4/AVALOKTMP/==MATLAB==/ActiveProjects/DNABIOMETRICS/PSUPITTMERGED/OUTPUTSNPJOBS/';
nrComp = 17;
for i=1:nrComp
   removeActiveFilesAfterShutDown('path',path,'id',i); 
end