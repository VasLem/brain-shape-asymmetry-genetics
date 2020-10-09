function out = applyScoreFusion(input,fuser,varargin)
         options = readVarargin(varargin{:});
         
         out.nF = input.nF;
         out.nS = 1;
         out.nGEN = input.nGEN;
         out.nIMPOS = input.nIMPOS;
         
         out.GEN = getPrediction(fuser,input.GEN',options.type)';
         out.IMPOS = getPrediction(fuser,input.IMPOS',options.type)';
         out.REDIMPOS = getPrediction(fuser,input.REDIMPOS',options.type)';
         
         if options.display % need to adapt to N-D
            
            figure;
            subplot(1,2,1);hist(out.GEN,30);title('GEN');
            subplot(1,2,2);hist(out.IMPOS,30);title('IMPOSTER');
            
            
            if input.nS>2, return;end
            if input.nS==1, return;end
            ngen = size(input.GEN,2);
            nimpos = size(input.REDIMPOS,2);
            class = [ones(1,ngen),-1*ones(1,nimpos)];
            FS = [input.GEN';input.REDIMPOS'];  
             
             
            h1 = (max(FS(:,1))-min(FS(:,1)))/100;
            h2 = (max(FS(:,2))-min(FS(:,2)))/100;
            [X1,X2] = meshgrid(min(FS(:,1)):h1:max(FS(:,1)),...
                               min(FS(:,2)):h2:max(FS(:,2)));
            [pred,score] = predict(fuser.Mdl,[X1(:),X2(:)]);
            binGrid = reshape(pred,size(X1,1),size(X2,2));
            
            val = score(:,2);
            scoreGrid = reshape(val,size(X1,1),size(X2,2));
            figure;sf = surf(X1,X2,zeros(size(X1)),scoreGrid,'LineStyle','none');view(0,90);
            hold on;
            ind = class==1;
            plot(FS(ind,1),FS(ind,2),'ko');
            ind = class==-1;
            plot(FS(ind,1),FS(ind,2),'ro');
             
         end
end

function options = readVarargin(varargin)
                 Input = find(strcmpi(varargin, 'options'));
                 if ~isempty(Input), options = varargin{Input+1}; return;end    
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),options.display = false;else, options.display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'type'));
                 if isempty(Input),options.type = 'soft';else, options.type = varargin{Input+1};end
end

function out = getPrediction(fuser,FS,type)
         [pred,score] = predict(fuser.Mdl,FS);
         switch type
             case 'hard'
                 out = pred;
             case 'soft'
                 out = score(:,2);
         end
end