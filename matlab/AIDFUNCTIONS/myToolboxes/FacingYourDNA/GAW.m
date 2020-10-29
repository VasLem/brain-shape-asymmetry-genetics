function [bestW,avgW] = GAW(Matches,wM,n,g,GAoptions)
         %GAoptions.PopSize = 100;
         % initialize
         nW = size(Matches,1);
         iter = 0;
         redind = ceil(GAoptions.Red*GAoptions.PopSize);
         % initialize
         Pop = rand(GAoptions.PopSize,nW);
         val = evalPop(Pop,Matches,wM,n,g);
         if GAoptions.Display, disp(['Iter: ' num2str(0) ' Best: ' num2str(min(val)) ' Mean: ' num2str(mean(val))]); end
         iterate = true;
         while iter<GAoptions.MaxIter&&iterate
            iter = iter + 1;
           % SORT
            [val,sortind] = sort(val); 
           % GENERATE
            Pop = nextPop(Pop,sortind,redind);
           % EVAL
            val(redind+1:end) = evalPop(Pop(redind+1:end,:),Matches,wM,n,g);
            if GAoptions.Display, disp(['Iter: ' num2str(iter) ' Best: ' num2str(min(val)) ' Mean: ' num2str(mean(val))]); end
            %if length(find(ranks==0))>= options.minRank0, iterate = false; end
         end
         bestV = min(val);
         bestW = mean(Pop(val==bestV,:),1)';
         if nargout < 2, return; end
         [~,ind] = sort(val);
         scaling = zeros(size(val));
         scaling(ind) = 1./ ((1:length(val)).^0.5);
         scaling = GAoptions.PopSize * scaling ./ sum(scaling);
         avgW = (sum(repmat(scaling',1,nW).*Pop)./sum(scaling))';
end

function val = evalPop(Pop,Matches,wM,n,g)
         val = nan*zeros(1,size(Pop,1));
         parfor i=1:1:length(val)
             [R,EER] = optWScoreAndEval(Matches,wM,Pop(i,:)',n,g);
             val(i) = R+EER;
         end
end

function Pop = nextPop(Pop,sortind,redind)
         Pop = Pop(sortind,:);PopSize = size(Pop,1);
         PopR = Pop(1:redind,:);
         RemainSize = PopSize - redind;
         h = 1.06*min(iqr(PopR)/1.349,std(PopR))./(redind^(1/5));
         PopR = PopR';c= size(PopR,2);
         randsel = bsxfun(@(x,y) x(randperm(y(1))),PopR',repmat(c,1,c)');
         Pop(redind+1:end,:) = randsel(1:RemainSize,:) + repmat(h,RemainSize,1).*randn(RemainSize,size(Pop,2));
         Pop(Pop<0) = 0;% exclude negative weigths
         %Pop(Pop>1) = 1;% upper boundary on weigths
end

