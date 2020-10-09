function a = normaliseArray(a)

    nrchannels = numel(size(a));
    if nrchannels == 3
       for i=1:1:size(a,3)
           mi = nanmin(nanmin(a(:,:,i)));
           ma = nanmax(nanmax(a(:,:,i)));
           a(:,:,i) = (a(:,:,i)-mi)./(ma-mi);
       end
    else
       mi = nanmin(a(:));
       ma = nanmax(a(:));
       a = (a-mi)./(ma-mi);
    end
end


