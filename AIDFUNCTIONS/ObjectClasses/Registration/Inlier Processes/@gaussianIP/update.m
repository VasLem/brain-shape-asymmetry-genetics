function out = update(obj,Tmodel)
% updating the standard deviation parameters
    if isempty(obj.Smeasure), error('Cannot update IP, no Smeasure set'); end    
    %eval(obj.Smeasure,Tmodel);
    [B,sumB] =  getB(obj,Tmodel);
    P = zeros(obj.nrSM,1);
    for i=1:1:obj.nrSM
        P(i) = sqrt(sum(B.*obj.Smeasure.Evaluation(i,:).^2)/sumB);
        P(i) = P(i) + obj.Temp*P(i);% Temp can increase STD width by 2 if Temp == 1;
    end
    if nargout == 1, out = P; return; end
    obj.PdfP = P;
end