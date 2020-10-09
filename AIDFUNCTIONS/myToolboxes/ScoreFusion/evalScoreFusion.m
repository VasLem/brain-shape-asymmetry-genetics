function out = evalScoreFusion(input,varargin)

out.VER.

[~,~,out.auc,X,Y,~,~,out.pauc] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));%#ok<*FNDSB> % CHECK ORDER
   if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--'); 
                    plot(X,Y,'b-','LineWidth',1.5);
                 end







end


function [EER,G,AUC,x,y,TH,Y,pAUC] = getEER(tmatches,fmatches)
                g = [ones(1,length(tmatches)), -1*ones(1,length(fmatches))];
                [sorted,order] = sort([tmatches;fmatches],'ascend');
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