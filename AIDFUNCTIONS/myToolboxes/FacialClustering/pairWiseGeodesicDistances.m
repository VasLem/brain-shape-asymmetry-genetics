function D = pairWiseGeodesicDistances(scan)
    D = zeros(scan.nrV,scan.nrV);
    [path,ID] = setupParForProgress(scan.nrV);
    parfor k=1:scan.nrV
        D(k,:) = intraDistances(scan,'VertexIndex',k);
        parfor_progress;
    end
    closeParForProgress(path,ID);
end