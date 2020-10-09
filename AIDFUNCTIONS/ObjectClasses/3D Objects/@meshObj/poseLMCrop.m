function out = poseLMCrop(obj,dist)
 if nargout == 1;
    obj = clone(obj);
    out = obj;
 end
 if isempty(obj.PoseLM), indicatePoseLM(obj); end
 N = KNN(obj,obj.PoseLM,1);
 distances = intraDistances(obj,'VertexIndex',N,'Type','edge');
 index = find(distances<dist);
 crop(obj,'VertexIndex',index);
end