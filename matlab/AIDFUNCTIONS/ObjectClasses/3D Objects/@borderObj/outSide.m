function out = outSide(obj,points,varargin)
         th = readVarargin(varargin{:});
         eval = rbfEval(obj,points,varargin);
         out = find(eval<th);
end

function th = readVarargin(varargin)
         Input = find(strcmp(varargin,'TH'));
         if ~isempty(Input)
            th = varargin{Input+1};
         else
            th = 0;
         end
end