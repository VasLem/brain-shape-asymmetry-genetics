function pvalue = myHWE(geno)
        tmp = geno(geno>=0);
        pvalue=hwetest([sum(tmp==0) sum(tmp==1) sum(tmp==2)]);  
end