function derivateCheck(obj,Tmodel,FS,Pindex)
            % initialize both
                initialize(Tmodel,FS);
                initialize(obj,Tmodel);
                eval(obj,Tmodel);
            % Derivate 
                Tmodel.ActiveP = Pindex;
                derivateField(Tmodel,FS,'All');
                derivate(obj,Tmodel);
                %D = obj.Derivative;           
            % initialize loop
                delta = Tmodel.d(Pindex);
                range = (-delta/2:delta/5:delta/2);
                eval(obj,Tmodel);
                values = zeros(1,length(range));
                grads = obj.Evaluation*ones(1,length(range));
                origP = Tmodel.P;
            % Do loop
                for i=1:1:length(range)
                    Tmodel.P = origP;
                    Tmodel.P(Pindex) = origP(Pindex) + delta*range(i);
                    eval(Tmodel,FS);
                    eval(obj,Tmodel);
                    values(i) = obj.Evaluation;
                    grads(i) = grads(i) + delta*range(i)*obj.Derivative;
                end
            % Reset Tmodel    
                Tmodel.P = origP;
            % Plot Result    
                figure; hold on;
                plot(values,'b-','LineWidth',1.5);
                plot(grads,'r-','LineWidth',1);            
end