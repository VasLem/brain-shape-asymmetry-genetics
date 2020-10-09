function out = extractGWASPeaks(GWAS,dist,LD,pCrit,display)
    if nargin<5, display = false; end     
    TEST = (GWAS.pValues<=pCrit)';
    TEST = sum(TEST,2);
    index = find(TEST);
    bestP = min(GWAS.pValues)';
    
    % PW LD computations
    SNPS = GWAS.SNP;
    SNPS = double(SNPS);
    SNPS(SNPS==-1) = nan;

    nHits = length(index);
    list = 1:nHits;
    labels = zeros(1,nHits);
    top = zeros(1,nHits);
    POS = GWAS.POS(index);
    chr = GWAS.CHR(index);
    listP = bestP(index);
    SNPS = SNPS(:,index);
    lab = 0;
    
    while ~isempty(list)&&lab<=length(index)
          lab = lab+1;
          nlist = length(list);
          [m,mind] = min(listP(list));
          bestind = list(mind);
          top(bestind) = 1;
          labels(bestind) = lab;
          list = setdiff(list,bestind);% take the best out;
          bestchr = chr(bestind);
          subindex = find(chr(list)==bestchr);
          if isempty(subindex), continue; end
          bestsnp = SNPS(:,bestind);
          subind = list(subindex);
          bestpos = POS(bestind);
          subpos = POS(subind);
          take = subpos<=(bestpos+dist)&subpos>=(bestpos-dist);
          if sum(take)==0, continue; end
          subind = subind(take);
          subsnp = SNPS(:,subind);
          CORR = nan*zeros(1,length(subind));
          indi = find(~isnan(bestsnp));
          for r=1:length(subind)
              inj = subsnp(:,r);
              indj = find(~isnan(inj));
              indr = intersect(indi,indj);
              CORR(r) = abs(corr(bestsnp(indr),inj(indr)));
          end
          subind = subind(CORR>=LD);
          labels(subind) = lab;
          list = setdiff(list,subind);          
    end
    nPeak = sum(top);
    peakindex = find(top);
    %figure;hist(labels,nPeak);
    peakrs = GWAS.RS(index(peakindex));
    peaksupport = zeros(1,nPeak);
    for i=1:1:nPeak
        peaksupport(i) = sum(labels==labels(peakindex(i))); 
    end
    %good = find(peaksupport>10);
    if display
       CHR = cell(1,23);
        for c=1:1:23
            ind = find(GWAS.CHR==c);
            CHR{c}.P = GWAS.pValues(:,ind)';
            CHR{c}.nrSNP = length(GWAS.RS(ind));
            CHR{c}.MAF = GWAS.Balances(1,ind);
            CHR{c}.RS = GWAS.RS(ind);
        end 
        f = myManhattanplotPaper(CHR,peakrs,1:63,[]);
    end
    
    out.RS = GWAS.RS(index);
    out.CHR = GWAS.CHR(index);
    out.POS = GWAS.POS(index);
    out.P = GWAS.pValues(:,index);
    out.MAF = GWAS.Balances(1,index);
    out.SNP = GWAS.SNP(:,index);
    out.A1 = GWAS.A1(index);
    out.A2 = GWAS.A2(index);
    out.LABELS = labels;
    out.TOP = top;
    out.nPeak = nPeak;
    out.PeakIndex = peakindex;
    out.PeakLabel = labels(peakindex);
    out.PeakSupport = peaksupport;
    out.PeakRS = peakrs;
    out.PeakCHR = out.CHR(peakindex);
    out.PeakPOS = out.POS(peakindex);
    out.PeakpValue = min(out.P(:,peakindex));
    out.pCrit = pCrit;
    out.Dist = dist;
    out.LD = LD;
   
end