function parforFileSaving(JOBID,COMPID,JOBSTATUS,in,path) %#ok<*INUSL>
        if nargin<5, path = [pwd '/'];end
        switch JOBSTATUS
            case 0
                % FILE TOUCHING: BLOCKING JOB
                parforTouchUnixToDisk([path num2str(JOBID) '_' num2str(COMPID) '_0.mat']);
                return;
            case {1 2}
                % FIRST SAVING
                save([path num2str(JOBID) '_' num2str(COMPID) '_0'],'in');
                % THEN RENAMING
                parforRenameFileUnix([path num2str(JOBID) '_' num2str(COMPID) '_0.mat'],[path num2str(JOBID) '_' num2str(COMPID) '_' num2str(JOBSTATUS) '.mat']);
        end
end