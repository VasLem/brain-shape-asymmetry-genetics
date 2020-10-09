function out = vectorCorr(v1,v2)
         v1 = v1/norm(v1);v2 = v2/norm(v2);
         T = v1'*v2;N = sqrt((v1'*v1)*(v2'*v2));
         out = T/N;
end