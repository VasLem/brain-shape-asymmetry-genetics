function out = eval(obj,Tmodel)
% negative log evaluation of gaussian pdf
    if isempty(obj.Smeasure), error('Cannot eval IP, no Smeasure set'); end
    if ~obj.FastEval
       eval(obj.Smeasure,Tmodel);
    end
    pdf(obj,Tmodel);
    [B,sumB] =  getB(obj,Tmodel); %#ok<NASGU>
    %warning('off','All');
    tmp = obj.PdfEvaluation;
    tmp = max(tmp,1e-100*ones(size(tmp)));% prevent Log(0) = -inf;
    %tmp = B.*log(obj.PdfEvaluation);
    %warning('on','All');
    %tmp(find(isinf(tmp))) = 10000; %#ok<FNDSB>
    %val = -1*sum(tmp);
    val = -1*sum(B.*log(tmp));
%    val = -1*sum(B.*log(obj.PdfEvaluation));
%     index = find(obj.PdfEvaluation);
%     val = -1*sum(B(index).*log(obj.PdfEvaluation(index)));
    if nargout == 1, out = val; return; end
    obj.Evaluation = val;
end