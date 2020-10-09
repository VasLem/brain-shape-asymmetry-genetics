function out = getAverageFromDirectory(obj,filenames)
    % out = getAverage(obj,in) or getAverage(obj,in)
    % Compute the Average for the PCA space based on give training data
    % INPUT
    % obj = PCA space object
    % in = Data matrix or Data structure
    % OUTPUT
    % out = new PCA space object or out = updated obj;
    %
    % created by Peter Claes
    if nargout == 1,obj = clone(obj);out = obj;end
    nrfiles = length(filenames);
    AvgVec = zeros(size(obj.AvgVec));
    f =waitbar(0,'Average f. files');drawnow;
    disp('Starting 2');
    counter = 0;
    for i=1:1:nrfiles
        try
            load(filenames{i});
            Vec = ass.Norm.TextureMap.Image(:)*255;
            if ~isempty(find(isnan(Vec))),Vec = ass.Scan.TextureMap.Image(:)*255;end  %#ok<*EFIND>
            if isempty(find(isnan(Vec)))
                AvgVec = AvgVec + Vec;
                counter = counter +1;
            end
            waitbar(i/nrfiles,f);drawnow;
        catch
            disp('Missing Texture');
        end
    end
    delete(f);
    obj.AvgVec = AvgVec/counter;
end