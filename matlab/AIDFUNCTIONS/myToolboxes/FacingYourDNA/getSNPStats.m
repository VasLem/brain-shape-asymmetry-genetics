function [MAF,CALL,M] = getSNPStats(in)
          tmp = double(in);
          N = size(tmp,1);
          Ne = N-sum(tmp==-1,1);
          CALL = Ne./N;
          tmp(tmp==-1)=0;
          Na = sum(tmp);
          p = Na./Ne./2;
          [MAF,ind] = min([p;1-p],[],1);
          switch ind
              case 1
                  M = 0;
              case 2
                  M = 2;
          end  
end