function matrix = CenterDataInMatrix(matrix,missing)
    if nargin<2, missing = false;end
    [~,N] = size(matrix);
    for i=1:N
        val = matrix(:,i);
        avg = nanmean(val);
        val = val-avg;
        if missing, val(isnan(val))=0; end
        matrix(:,i) = val;
    end
end