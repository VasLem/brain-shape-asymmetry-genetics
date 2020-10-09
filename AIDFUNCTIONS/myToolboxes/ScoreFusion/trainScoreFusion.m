function out = trainScoreFusion(input,type,varargin)
         options = readVarargin(varargin{:});
         out.Type = type;
         
         ngen = size(input.GEN,2);
         
%          nimpos = size(input.REDIMPOS,2);
%          class = [ones(1,ngen),-1*ones(1,nimpos)];
%          FS = [input.GEN';input.REDIMPOS'];      
         nimpos = size(input.IMPOS,2);
         class = [ones(1,ngen),-1*ones(1,nimpos)];
         FS = [input.GEN';input.IMPOS'];
         %Wgen = 1-ngen/(ngen+nimpos);
         %Wimpos = 1-nimpos/(ngen+nimpos);
         %W = [Wgen*ones(ngen,1);Wimpos*ones(nimpos,1)];
         W = [(1/ngen)*ones(ngen,1);(1/nimpos)*ones(nimpos,1)]/2;
         
         [i,j] = find(isnan(FS));
         index = setdiff(1:size(FS,1),unique(i));
         
         FS = FS(index,:);
         class = class(index);
         W = W(index);
         
         %W = ones(length(class),1);
         
         switch type
             case 'LDA'
                 Mdl = fitcdiscr(FS,class,'Weights',W,'DiscrimType','linear');
             case 'SVM'
                 Mdl = fitcsvm(FS,class,'KernelFunction','rbf','Weights',W,'KernelScale','auto');
             case 'TREE'
                 Mdl = fitctree(FS,class,'Weights',W);
             case 'KNN'
                 Mdl = fitcknn(FS,class,'Weights',W,'NumNeighbors',options.NumNeighbors,'Distance','euclidean');
             case 'NB'
                 Mdl = fitcnb(FS,class,'Weights',W);
         end
         out.Mdl = Mdl;
         out.W = W;
         if options.display
            if size(FS,2)>2, return;end
            h1 = (max(FS(:,1))-min(FS(:,1)))/100;
            h2 = (max(FS(:,2))-min(FS(:,2)))/100;
            [X1,X2] = meshgrid(min(FS(:,1)):h1:max(FS(:,1)),...
                               min(FS(:,2)):h2:max(FS(:,2)));
            [pred,score] = predict(Mdl,[X1(:),X2(:)]);
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
                 Input = find(strcmpi(varargin, 'test'));
                 if isempty(Input),options.test = [];else, options.test = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'NumNeighbors'));
                 if isempty(Input),options.NumNeighbors = 20;else, options.NumNeighbors = varargin{Input+1};end
end