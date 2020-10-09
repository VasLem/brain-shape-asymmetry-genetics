function out = getRVLMv2(RF)
     [nLM,nV,~] = size(RF);
     RV = zeros(nLM,nLM);
     disp('BUILDING RV MATRIX');
     [path,ID] = setupParForProgress(nLM);
     parfor i=1:nLM
        % i=1;
        datai = squeeze(RF(i,:,:))';
        tmp = zeros(1,nLM);
        for j=i:nLM
            %j=2;
            dataj = squeeze(RF(j,:,:))';
            covMatrix = cov([datai,dataj]);
            numRV = trace(covMatrix(1:nV,nV+1:end)*covMatrix(nV+1:end,1:nV)); % ok
            var1sum = sum(sum(covMatrix(1:nV,1:nV).^2));
            var2sum = sum(sum(covMatrix(nV+1:end,nV+1:end).^2));
            denRV = sqrt(var1sum*var2sum);
            tmp(j) = numRV/denRV;
        end
        RV(i,:) = tmp;
        parfor_progress;
     end
     closeParForProgress(path,ID);
     out = triu(RV,1)'+ triu(RV);
end