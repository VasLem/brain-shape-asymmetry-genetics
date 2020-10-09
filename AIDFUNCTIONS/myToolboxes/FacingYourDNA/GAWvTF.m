function [bestW,avgW] = GAWvTF(TMatches,FMatches,wM,GAoptions)
         %GAoptions.PopSize = 100;
         % initialize
         nW = size(wM,1);
         iter = 0;
         redind = ceil(GAoptions.Red*GAoptions.PopSize);
         % initialize
         Pop = rand(GAoptions.PopSize,nW);
         val = evalPop(Pop,TMatches,FMatches,wM);
         if GAoptions.Display, disp(['Iter: ' num2str(0) ' Best: ' num2str(min(val)) ' Mean: ' num2str(mean(val))]); end
         iterate = true;
         while iter<GAoptions.MaxIter&&iterate
            iter = iter + 1;
           % SORT
            [val,sortind] = sort(val); 
           % GENERATE
            Pop = nextPop(Pop,sortind,redind);
           % EVAL
            val(redind+1:end) = evalPop(Pop(redind+1:end,:),TMatches,FMatches,wM);
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

function val = evalPop(Pop,TMatches,FMatches,wM)
         val = nan*zeros(1,size(Pop,1));
         parfor i=1:1:length(val)
             [R,EER] = myWScoreAndEvalTF(TMatches,FMatches,wM,Pop(i,:)');
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

function [R,EER,ST,SF] = myWScoreAndEvalTF(TMatches,FMatches,wM,W)
         nF = size(FMatches,3);
         nT = size(TMatches,2);
         tmpW = repmat(W,1,nT).*wM;
         ST = getScores(TMatches,tmpW)';
         SF = getScores(FMatches,repmat(tmpW,1,1,nF));
         [R,EER] = evalScores(ST,SF,nT,nF); 
         %R = 0;EER = 0;
end

function Score = getScores(Matches,tmpW)
         Score = squeeze(nansum(tmpW.*Matches,1)./nansum(tmpW,1));
end

function [R,EER] = evalScores(ST,SF,nT,nF)
         % Rank analysis
         R = mean((sum(SF<repmat(ST,1,nF),2)+1)/(nF+1));
         % EER analysis
         %if nargout < 2, return; end
         g = [ones(1,nT), zeros(1,nT*nF)];
         [~,order] = sort([ST;SF(:)]);
         g = g(order);
         true_neg = g == 0;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
         true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
         y = tpf(1:end-1)+dy./2;
         x = fpf(1:end-1)+dx./2;
         yn = 1-y;d = abs(x-yn);
         [~,indmin] = min(d);
         EER = ((x(indmin)+yn(indmin))/2);
end

