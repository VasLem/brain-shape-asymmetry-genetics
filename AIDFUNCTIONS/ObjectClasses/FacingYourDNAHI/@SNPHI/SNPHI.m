classdef SNPHI < PPMHI
   % General Properties
   properties
      VarName = [];
      TestInd = 4;
      CHR = 0;
      POS = 0;
   end
   properties (Hidden = true)
      NanVal = -1;
      ALLELE1 = 'A1';
      ALLELE2 = 'A2';
   end
   properties (Dependent = true)
      RS;
      MAF;
      GT;
      GTFREQ;
      GTPROB;
      GTN;
      MinorityClass;
      TestLabel;
      mALLELE;
      MALLELE;
   end
   properties (Hidden = true, Dependent = true)
      VarType;
      ID;
      codeAL;
      nAL;
   end
   methods % CONSTRUCTOR
        function obj = SNPHI(varargin)
            obj = obj@PPMHI(varargin{:});         
        end
   end
   methods % GETTING/SETTING
       function out = get.RS(obj)
          out = cleanUpRS(obj.VarName); 
       end
       function set.RS(obj,in)
          obj.VarName = in; 
       end
       function out = get.GT(obj)
           out = obj.XX;
       end
       function set.GT(obj,in)
           obj.XX = in;
       end 
       function out = get.TestLabel(obj)
           switch obj.TestInd
               case 0
                   out = 'NO EFFECT';
               case 1
                   %out = 'RECESSIVE';
                   out = 'aa-aAAA';
               case 2
                   out = 'aaaA-AA';
               case 3
                   out = 'aaAA-aA';
               case 4
                   out = 'aa-aA-AA';
               case 5
                   out = 'aa-AA';
               case 6
                   out = 'aa-aA';
               case 7
                   out = 'aA-AA';
               case 8
                   out = 'NON ADDITIVE';
               case 9
                   out = 'ADDITIVE AA AB';
               case 10
                   out = 'ADDITIEVE AB BB';    
           end
       end
       function out = get.ID(obj) %#ok<*MANU>
           out = 1000;
       end
       function out = get.VarType(obj)
          out = 'SNP'; 
       end
       function out = get.MinorityClass(obj)
                out = 1;
       end
       function out = get.MAF(obj)
           tmp = obj.XG;
           out = (2*sum(tmp==1)+sum(tmp==0))/obj.nAL;
       end
       function out = get.codeAL(obj)
           if isempty(obj.XX), out = 2;end
           switch sum(obj.XX)>obj.nAL/2
               case true
                   out = 'A';
               case false
                   out = 'a';
           end
       end
       function out = get.nAL(obj)
          out = obj.nSNotNan*2; 
       end
       function out = get.GTFREQ(obj)
                tmp = obj.XG;
                out = [sum(tmp==1) sum(tmp==0) sum(tmp==-1)]./obj.nSNotNan;
       end
       function out = get.GTN(obj)
                tmp = obj.XG;
                out = [sum(tmp==1) sum(tmp==0) sum(tmp==-1)];
       end
       function out = get.GTPROB(obj)
                if isempty(obj.XX),out= [];return;end 
                pm = obj.MAF;
                pM = (1-pm);
                out = [pm*pm pm*pM pM*pM]; 
       end
       function out = get.mALLELE(obj)
           switch obj.codeAL
               case 'A'
                   out = obj.ALLELE1;
               case 'a'
                   out = obj.ALLELE2;
           end
       end
       function out = get.MALLELE(obj)
           switch obj.codeAL
               case 'a'
                   out = obj.ALLELE1;
               case 'A'
                   out = obj.ALLELE2;
           end
       end
   end
   methods % DATA CHECKING
   end
   methods % INTERFACING
       function out = GT2GG(obj,in)
           out = XG2GG(obj,XX2XG(obj,in));
       end
       function out = XX2XG(obj,in)
                 out = single(in);
                 switch obj.codeAL
                     case 'A'
                         out(in==2) = -1;
                         out(in==0) = 1;
                     case 'a'
                         out(in==2) = 1;
                         out(in==0) = -1;
                 end
                 out(in==1) = 0;
                 out(in==-1) = nan;
        end
       function out = XG2GG(obj,in) %#ok<*INUSL>
                 out = in;% Initialize
                 mmind = find(in==1); %#ok<*PROPLC>
                 Mmind = find(in==0);
                 MMind = find(in==-1);
                 switch obj.TestInd
                       case 0 % NO EFFECT
                           out = nan*out;
                       case 1 % test 1: AAAB BB (RECESSIVE)
                           out(union(Mmind,MMind)) = -1;
                           out(mmind) = 1;
                       case 2 % test 2: AA ABBB (DOMINANT)
                           out(MMind) = -1;
                           out(union(mmind,Mmind)) = 1;
                       case 3 % test 3: AABB AB (OVER DOMINANT)
                           out(union(mmind,MMind)) = -1;
                           out(Mmind) = 1;
                       case 4 % test 4 AA AB BB (ADDITIVE)
                           out(mmind) = 1;
                           out(Mmind) = 0;
                           out(MMind) = -1;
                       case 5 % test 5: AA BB
                           out(mmind) = 1;
                           out(Mmind) = nan;
                           out(MMind) = -1;
                       case 6 % test 6: AA AB
                           out(mmind) = 1;
                           out(Mmind) = -1;
                           out(MMind) = nan;
                       case 7 % test 7: AB BB
                           out(mmind) = nan;
                           out(Mmind) = 1;
                           out(MMind) = -1;
                end
        end
   end
   methods % CLASSIFIER
       function [G,COST,W,ind,Tind,Find] = prepTrainFITC(obj,X) %#ok<*INUSL>
            G = XX2GG(obj,X);
            ind = find(~(X==obj.NanVal));
            Find = find(G==-1*obj.MinorityClass);
            Tind = find(G==obj.MinorityClass);
            COST.ClassNames = [obj.MinorityClass -1*obj.MinorityClass];
            COST.ClassificationCosts = zeros(2,2);
            COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
            COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
            W = nan*zeros(size(X));
            W(Tind) = 1-(length(Tind)/length(ind));
            W(Find) = 1-(length(Find)/length(ind));
       end
       function [G,COST,W,ind,Tind,Find] = prepTrainSVM(obj,X) %#ok<*INUSL>
            G = XX2GG(obj,X);
            ind = find(~(X==obj.NanVal));
            Find = find(G==-1*obj.MinorityClass);
            Tind = find(G==obj.MinorityClass);
            COST.ClassNames = [obj.MinorityClass -1*obj.MinorityClass];
            COST.ClassificationCosts = zeros(2,2);
            COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
            COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
            W = nan*zeros(size(X));
            W(Tind) = 1-(length(Tind)/length(ind));
            W(Find) = 1-(length(Find)/length(ind));
       end
       function out = getCLASSGTW(obj,X)
           %options = PPMHI.readVarargin(varargin{:});
                 XG = XX2XG(obj,X);
                 GG = XX2GG(obj,X);
                 n = length(X);
                 out = nan*zeros(1,n,n);
                 forout = nan*zeros(n,n);
                 for j=1:1:n
                     switch XG(j)
                         case -1
                             forout(:,j) = 1-obj.GTFREQ(1);
                         case 0
                             forout(:,j) = 1-obj.GTFREQ(2);
                         case 1
                             forout(:,j) = 1-obj.GTFREQ(3);
                         otherwise
                     end
                     %index = find(obj.CLASS{ind}.MDL.ClassNames==G(j));
                     %if ~isempty(index), forout(:,j) = posterior(:,index); end
                 end
                 out(1,:,:) = forout;
       end
   end
   methods (Static = true)
   end
end