function parforAlive(COMPID)
         if COMPID==1
            parfor_progress;
         else
            time = clock;
            disp([date ' Hour: ' num2str(time(4)) ':' num2str(time(5)) ':' num2str(round(time(6)))]);
         end
end