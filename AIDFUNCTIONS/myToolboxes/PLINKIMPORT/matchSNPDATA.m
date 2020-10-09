function [out] = matchSNPDATA(in1,in2,type)
         if nargin<3, type = 'intersect';end
         
         list1 = [string(in1.POS), string(in1.A1), string(in1.A2)];   
         list2 = [string(in2.POS), string(in2.A1), string(in2.A2)];
         list2F = [string(in2.POS), string(in2.A2), string(in2.A1)];
             
         [~,ind12] = ismember(list1,list2,'rows');
         ind21 = find(ind12);
         ind12 = ind12(ind21);
         nO = length(ind12);
             
         [~,ind12F] = ismember(list1,list2F,'rows');
         ind21F = find(ind12F);
         ind12F = ind12F(ind21F);
         [ind21F,IA] = setdiff(ind21F,ind21);
         ind12F = ind12F(IA);
         nF = length(ind12F);               
         
         
         
         n = nO + nF;
                     SNP = [[SNP1(:,ind21); SNP2(:,ind12)],[SNP1(:,ind21F); SNP2F(:,ind12F)]];
                     POS = [in1.POS(chrind1(ind21));in1.POS(chrind1(ind21F))];
                     A1 = [in1.A1(chrind1(ind21));in1.A1(chrind1(ind21F))];
                     A2 = [in1.A2(chrind1(ind21));in1.A2(chrind1(ind21F))];
                     RSID = [in1.RSID(chrind1(ind21));in1.RSID(chrind1(ind21F))];
                     [POS,sortindex] = sort(POS,'ascend');
                     SNP = SNP(:,sortindex);
                     nSNP = size(SNP,2);
                     A1 = A1(sortindex);
                     A2 = A2(sortindex);
                     RSID = RSID(sortindex);  
         
         
         for c=1:1:nchrs
             % c = 1;
             chrind1 = find(in1.CHR==chrs(c));
             chrind2 = find(in2.CHR==chrs(c));

             list1 = [string(in1.POS(chrind1)), string(in1.A1(chrind1)), string(in1.A2(chrind1))];
             
             list2 = [string(in2.POS(chrind2)), string(in2.A1(chrind2)), string(in2.A2(chrind2))];
             list2F = [string(in2.POS(chrind2)), string(in2.A2(chrind2)), string(in2.A1(chrind2))];
             
             [~,ind12] = ismember(list1,list2,'rows');
             ind21 = find(ind12);
             ind12 = ind12(ind21);
             nO = length(ind12);
             
             [~,ind12F] = ismember(list1,list2F,'rows');
             ind21F = find(ind12F);
             ind12F = ind12F(ind21F);
             [ind21F,IA] = setdiff(ind21F,ind21);
             ind12F = ind12F(IA);
             nF = length(ind12F);                                  
%              
%              [test,IA,IB] = intersect(ind21,ind21F);
%              
%              
%              in1.RSID(ind21F(IB(1)))
%              in1.POS(ind21F(IB(1)))
%              in1.A1(ind21F(IB(1)))
%              in1.A2(ind21F(IB(1)))
%              in1.MAF(ind21F(IB(1)))
%              
%              in2.RSID(ind12F(IB(1)))
%              in2.POS(ind12F(IB(1)))
%              in2.A1(ind12F(IB(1)))
%              in2.A2(ind12F(IB(1)))
%              in2.MAF(ind12F(IB(1)))
%              
%              
%              in1.RSID(ind21(IA(1)))
%              in1.POS(ind21(IA(1)))
%              in1.A1(ind21(IA(1)))
%              in1.A2(ind21(IA(1)))
%              
%              in2.RSID(ind12(IA(1)))
%              in2.POS(ind12(IA(1)))
%              in2.A1(ind12(IA(1)))
%              in2.A2(ind12(IA(1)))
%              in2.MAF(ind12(IA(1)))                  
             
             SNP1 = in1.SNP(:,chrind1); 
             SNP2 = in2.SNP(:,chrind2);
             SNP2F = SNP2;SNP2F(SNP2==2) = 0;SNP2F(SNP2==0) = 2;
             switch type
                 case 'intersect'
                     n = nO + nF;
                     SNP = [[SNP1(:,ind21); SNP2(:,ind12)],[SNP1(:,ind21F); SNP2F(:,ind12F)]];
                     POS = [in1.POS(chrind1(ind21));in1.POS(chrind1(ind21F))];
                     A1 = [in1.A1(chrind1(ind21));in1.A1(chrind1(ind21F))];
                     A2 = [in1.A2(chrind1(ind21));in1.A2(chrind1(ind21F))];
                     RSID = [in1.RSID(chrind1(ind21));in1.RSID(chrind1(ind21F))];
                     [POS,sortindex] = sort(POS,'ascend');
                     SNP = SNP(:,sortindex);
                     nSNP = size(SNP,2);
                     A1 = A1(sortindex);
                     A2 = A2(sortindex);
                     RSID = RSID(sortindex);  
                 case 'add'
                     nSNP = size(SNP1,2);
                     tmpSNP = -1*ones(nS2,nSNP,'int8');
                     tmpSNP(:,ind21) = SNP2(:,ind12);
                     tmpSNP(:,ind21F) = SNP2F(:,ind12F);
                     SNP = [SNP1;tmpSNP];
                     POS = in1.POS(chrind1);
                     A1 = in1.A1(chrind1);
                     A2 = in1.A2(chrind1);
                     RSID = in1.RSID(chrind1);
             end
             [MAF,CALL,~,~,~] = getSNPStat(SNP);
             out.CHR = [out.CHR c*ones(nSNP,1)];
             out.RSID = [out.RSID; RSID(:)];
             out.POS = single([out.POS; POS(:)]);
             out.SNP = [out.SNP, SNP];
             out.CALL = [out.CALL; CALL(:)];
             out.MAF = [out.MAF; MAF(:)];
             out.A1 = [out.A1; A1(:)];
             out.A2 = [out.A2; A2(:)];
             out.nSNP = out.nSNP+nSNP;
         end
end

function [MAF,CALL,N,Ne,p] = getSNPStat(Genotypes)
         N = size(Genotypes,1);
         Ne = N-sum(Genotypes==-1,1);
         CALL = Ne./N;
         tmp = Genotypes;
         tmp(Genotypes==-1)=0;
         Na = sum(tmp);
         p = Na./Ne./2;
         MAF = min([p;1-p],[],1);
end