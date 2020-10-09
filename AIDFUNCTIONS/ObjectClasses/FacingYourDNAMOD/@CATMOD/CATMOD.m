classdef CATMOD < PPMMOD
    properties
    end
    properties (Dependent = true)
    end
    properties (Hidden = true) 
    end
    properties(Hidden = true, Dependent = true)
    end
    methods % CONSTRUCTOR
        function obj = CATMOD(varargin)
            obj = obj@PPMMOD(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING     
    end
    methods % BIOMETRICS
        function [EER,G,AUC,x,y,R,PREC,REC,GH,TNF,TH,Y,XH,YH] = oppBiometricTest(obj,TrX,TrRIP,TestX,TestRIP,posclass)
                n = length(TestRIP);
                C = unique(TrX);nrC = length(C);allC = 1:nrC;
                AVG = nan*zeros(1,nrC);STD = nan*zeros(1,nrC);
                for c=1:1:nrC
                    trindex = find(TrX==C(c));
                    %AVG(c) = mean(TrRIP(trindex));
                    %STD(c) = std(TrRIP(trindex));
                    [AVG(c),STD(c),MED,MAD] = PPMMOD.getRIPStatistics(TrRIP(trindex)',obj.Kappa);
                    %[~,~,AVG(c),STD(c)] = PPMMOD.getRIPStatistics(TrRIP(trindex)',obj.Kappa);
                end
                PDF = normpdf(repmat(TestRIP,1,nrC),repmat(AVG,n,1),repmat(STD,n,1));
                % true matches only defined by minority class
                cp = find(C==posclass);
                indexp = find(TestX==posclass);
                indexn = setdiff(1:n,indexp);
                remCp = setdiff(allC,cp);
                switch obj.Match
                    case 'QD'
                        bel = (PDF(:,cp))./PDF(:,remCp);
                        ratio = -log(bel);
                        %th = [min(bel):0.01:max(bel)];
                        th = 1;
                    case 'INLIER'
                        bel = (PDF(:,cp))./sum(PDF,2);
                        ratio = 1-bel;
                        %th = [0:0.01:1];
                        th = 0.5;
                end
                [EER,G,AUC,x,y,TH,Y] = PPMMOD.getEER(ratio(indexp),ratio(indexn));
                %figure;plot(x,y,'b-');
                R = PPMMOD.getRANK(ratio(indexp),ratio(indexn));
                PREC = nan*zeros(1,length(th));
                REC = nan*zeros(1,length(th));
                GH = nan*zeros(1,length(th));
                XH = nan*zeros(1,length(th));
                YH = nan*zeros(1,length(th));
                %GH(2,:) = th;
                TNF = nan*zeros(1,length(th));
                for i=1:1:length(th)
                    Tclass = bel(indexp)>=th(i);
                    Fclass = bel(indexn)>=th(i);
                    [PREC(i),REC(i),GH(1,i),TNF(i),XH(i),YH(i)] = PPMMOD.getHardClass(Tclass,Fclass);
                end
        end
    end
    methods % INTERFACING
    end
    methods (Static = true)    
    end
end


% function [EER,G,AUC,R] = oppBiometricTest(obj,TrX,TrRIP,TestX,TestRIP,posclass)
%                 n = length(TestRIP);
%                 C = unique(TrX);nrC = length(C);allC = 1:nrC;
%                 AVG = nan*zeros(1,nrC);STD = nan*zeros(1,nrC);
%                 for c=1:1:nrC
%                     trindex = find(TrX==C(c));
%                     %AVG(c) = mean(TrRIP(trindex));
%                     %STD(c) = std(TrRIP(trindex));
%                     [AVG(c),STD(c),MED,MAD] = PPMMOD.getRIPStatistics(TrRIP(trindex)',obj.Kappa);
%                     %[~,~,AVG(c),STD(c)] = PPMMOD.getRIPStatistics(TrRIP(trindex)',obj.Kappa);
%                 end
%                 PDF = normpdf(repmat(TestRIP,1,nrC),repmat(AVG,n,1),repmat(STD,n,1));
%                 % true matches only defined by minority class
%                 cp = find(C==posclass);
%                 indexp = find(TestX==posclass);
%                 indexn = setdiff(1:n,indexp);
%                 remCp = setdiff(allC,cp);
%                 switch obj.Match
%                     case 'QD'
%                         ratio = (PDF(:,cp))./PDF(:,remCp);
%                         switch obj.Classify
%                             case 'SOFT'
%                                 ratio = -log(ratio);
%                             case 'HARD'
%                                 ratio = ratio<1;
%                         end
%                     case 'INLIER'
%                         ratio = (PDF(:,cp))./sum(PDF,2);
%                         switch obj.Classify
%                             case 'SOFT'
%                                 ratio = 1-ratio;
%                             case 'HARD'
%                                 ratio = ratio<0.5;
%                         end
%                 end
%                 [EER,G,AUC,x,y] = PPMMOD.getEER(ratio(indexp),ratio(indexn));
%                 figure;plot(x,y,'b-');
%                 R = PPMMOD.getRANK(ratio(indexp),ratio(indexn));
% end