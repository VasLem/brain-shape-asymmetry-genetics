function out = derivate(obj,Tmodel)
         % TM needs to be evaluated before SM derivate !!!!
         %if isempty(Tmodel.Evaluation),error('Tmodel not evaluated');end
         if ~isempty(Tmodel.VertDerivative)
             %disp('Existing vertD');
             % SM needs to be evaluated before SM derivative
             % only simple test, it is adviced to always eval once before derivates in Main
             % procedure!!!!
             %if ~(size(obj.rbfGrad,2)==Tmodel.nrV), eval(obj,Tmodel); end
             dist = max(obj.Distance,1e-100*ones(size(obj.Distance)));% robust against dividing by zero
             Pgrad = (1./dist).*sum(obj.Difference.*Tmodel.VertDerivative);
         else
             Pgrad = zeros(1,Tmodel.nrV);
         end
         if ~isempty(Tmodel.DistDerivative)
             %disp('Existing distD');
             Pgrad = Pgrad-Tmodel.DistDerivative;
         end       
         if nargout == 1,out = Pgrad; return; end
         obj.Derivative = Pgrad;
end