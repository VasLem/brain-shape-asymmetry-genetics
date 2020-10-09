function parforFileSavingv2(ID,COMPID,JOBSTATUS,nH,in) %#ok<INUSD>
        switch JOBSTATUS
            case 0
                % FILE TOUCHING: BLOCKING JOB
                parforTouchUnixToDisk([ID '_' num2str(COMPID) '_0.mat']);
                return;
            case {1 2}
                % FIRST SAVING
                save([ID '_' num2str(COMPID) '_0.mat'],'in');
                % THEN RENAMING
                parforRenameFileUnix([ID '_' num2str(COMPID) '_0.mat'],[ID '_' num2str(nH) '_' num2str(COMPID) '_' num2str(JOBSTATUS) '.mat']);
                %movefile([num2str(JOBID) '_' num2str(COMPID) '_0.mat'],[num2str(JOBID) '_' num2str(COMPID) '_' num2str(JOBSTATUS) '.mat']); 
        end
end