function out = evalBiometricMatch(match,varargin)
         nF = size(match,2);
         nS = size(match,1);
         options = readVarargin(varargin{:});
         if (options.standardize&&nS>1)
            for s=1:1:nS
               % do something
               tmp = squeeze(match(s,:,:));
               mintmp = min(tmp(:));
               tmp = tmp-mintmp;
               maxtmp = max(tmp(:));
               tmp = tmp./maxtmp;
               match(s,:,:) = tmp;
            end
         end
         if nS>1
            match = nanmean(match,1);% simple mean rule     
         end
         match = squeeze(match);
         ID = eye(nF,nF);
         indgen = find(ID==1);ngen = length(indgen);out.nGEN = ngen;
         indimpos = find(ID==0);nimpos = length(indimpos);out.nIMPOS = nimpos;
         % Verficiation analysis
         [out.VER.EER,out.VER.G,out.VER.auc,out.VER.x,out.VER.y,out.VER.TH,out.VER.Y,out.VER.pauc] = ...
                                                            getEER(match(indgen),match(indimpos));
         [out.IDEN.R,out.IDEN.CUMR] = getRANK(match);
         out.LINEUP = lineUps(match);
         if options.display
            if isempty(options.idenfigure) 
               options.idenfigure = figure;
               hold on;title('Identification');
               set(gca,'ylim',[0 100],'xlim',[1 100]);
               xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
               plot(1:100,1:100,'k--');
            end
            figure(options.idenfigure);
            plot(1:100,out.IDEN.CUMR,'Color',options.color,'LineWidth',1.5); 
            if isempty(options.verfigure)
               options.verfigure = figure;
               hold on;title('Verification');
               xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
               plot(0:0.01:1,0:0.01:1,'k--');
               plot(0:0.01:1,1:-0.01:0,'k-');
            end
            figure(options.verfigure);
            plot(out.VER.x,out.VER.y,'Color',options.color,'LineWidth',1.5);
         end
         out.options = options;
end

function options = readVarargin(varargin)
                 Input = find(strcmpi(varargin, 'options'));
                 if ~isempty(Input), options = varargin{Input+1}; return;end    
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),options.display = false;else, options.display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'standardize'));
                 if isempty(Input),options.standardize = false;else, options.standardize = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'verfigure'));
                 if isempty(Input),options.verfigure = [];else, options.verfigure = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'idenfigure'));
                 if isempty(Input),options.idenfigure = [];else, options.idenfigure = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'color'));
                 if isempty(Input),options.color = 'b-';else, options.color = varargin{Input+1};end
end


function [EER,G,AUC,x,y,TH,Y,pAUC] = getEER(tmatches,fmatches)
                g = [ones(1,length(tmatches)), -1*ones(1,length(fmatches))];
                [sorted,order] = sort([tmatches;fmatches],'descend');
                g = g(order);
                true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                y = tpf(1:end-1)+dy./2;
                x = fpf(1:end-1)+dx./2;
                [Y,indy] = max(y-x);
                yn = 1-y;d = abs(x-yn);
                [~,ind] = min(d);
                EER = ((x(ind)+yn(ind))/2);
                AUC = sum(dx.*y);
                FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
                TP = na-FN;FP = nn-TN;
                G = 1-sqrt((TP/(TP+FN))*(TN/(TN+FP)));
                T(1,1) = TP;
                T(1,2) = FP;
                T(2,1) = FN;
                T(2,2) = TN;
                TH(1) = 1-sorted(ind);
                TH(2) = 1-sorted(indy);
                se = PPMMV2.standardError(AUC,na,nn);
                if AUC>=0.5
                   pAUC = normpdf((AUC-0.5)/se,0,1);
                else
                   pAUC = 1;
                end
end

function [R,cumrank] = getRANK(match)
         nF = size(match,1);
         R = nan*zeros(1,nF);
         for f=1:1:nF
             %indrest = setdiff(1:nF,f);
             R(f) = sum(match(:,f)>=match(f,f)); 
         end
         R = (R./nF)*100;
         cumrank = zeros(1,100);
         for cr=1:1:100
             cumrank(cr) = (sum(R<=cr)./nF)*100;
         end
end

function out = lineUps(match)
         poolsizes = [2 5 10 20 40 100];
         npool = length(poolsizes);
         out = zeros(10,npool);
         nS = size(match,1);
         for r=1:1:10
             for n=1:npool
                 ps = poolsizes(n)-1;
                 succ = zeros(1,nS);
                 for t=1:nS
                     ind = randsample(setdiff(1:nS,t),ps);
                     fmatches = match(ind,t);
                     if sum(fmatches>=match(t,t))==0, succ(t) = 1; end
                 end
                 out(r,n) = (sum(succ)/nS)*100;
             end
         end
         out = mean(out);
end


