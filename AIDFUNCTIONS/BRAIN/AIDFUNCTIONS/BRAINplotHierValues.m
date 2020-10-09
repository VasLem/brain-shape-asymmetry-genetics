function out = BRAINplotHierValues(rend,SegmentVal,nLevels,varargin)
    values = zeros(1,rend.UHI.nLC);
    values(find(rend.UMASK)) = SegmentVal;
    out = plotHierValues3DwMASK(values,rend.UMASK,'nLevels',nLevels,varargin{:});
end