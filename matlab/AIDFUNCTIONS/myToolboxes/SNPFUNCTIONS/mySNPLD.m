function out = mySNPLD(geno1,geno2,model)
% THIS IS EM BASED LD COMPUTATION FOLLOWING THE PAPER:
% Hill (1974), Estimation of Linkage Disequilibrium in randomly mating
% populations
         % input
         if nargin<3, model = 'codominant-codominant'; end
         % define overlapping samples
         index = intersect(find(~isnan(geno1)),find(~isnan(geno2)));
         geno1 = geno1(index);geno2 = geno2(index);     
         % code towards the common allele
         geno1 = codetoM(geno1);
         geno2 = codetoM(geno2);
         switch lower(model)
             case 'codominant-codominant'
                 % codominant-codominant
                 Table = rot90(crosstab(geno1,geno2),2);
                 Nrow = zeros(1,3);
                 Ncol = zeros(1,3);
                 for i=1:1:3
                     Nrow(i) = sum(Table(i,:));
                     Ncol(i) = sum(Table(:,i));
                 end
                 N = sum(Nrow);
                 % Derived totals
                 X11 = 2*Table(1,1)+Table(1,2)+Table(2,1);
                 X21 = 2*Table(3,1)+Table(2,1)+Table(3,1);
                 X12 = 2*Table(1,3)+Table(1,2)+Table(2,3);
                 X22 = 2*Table(3,3)+Table(2,3)+Table(3,2);
                 p = (Nrow(1)+0.5*Nrow(2))/N;
                 q = (Ncol(1)+0.5*Ncol(2))/N; 
                 f11 = ((4*N)^-1)*(X11-X12-X21+X22)+0.5-((1-p)*(1-q));
                 acc = 1;
                 counter = 1;
                 while acc>1e-8&&counter<100
                    f11old = f11;
                    f11 = (X11+(Table(2,2)*f11*(1-p-q+f11))/(f11*(1-p-q+f11)+(p-f11)*(q-f11)))/(2*N);     
                    acc = abs(f11-f11old);          
                    counter = counter+1; 
                 end
                 D = f11-p*q;
                 switch sign(D)
                     case -1
                         Dmax = max([-1*p*q -1*(1-p)*(1-q)]);
                     case 1
                         Dmax = min([p*(1-q) (1-p)*q]);
                 end
                 r = D/sqrt(p*(1-p)*q*(1-q));
                 k = (N*D^2)/(p*(2-p)*q*(2-q));
                 R2 = k/N;
                 pval = 1-chi2cdf(k,1);
                 out.D = D;
                 out.Dmax = Dmax;
                 out.Dprime = D/Dmax;
                 out.k = k;
                 out.R2 = R2;
                 out.r = r;
                 out.r2 = r^2;
                 out.p = pval;
             case 'dominant-codominant'
                 % more complicated
                 out = [];
             case 'dominant-dominant'
                 % more complicated
                 out = [];
             case 'all'
                 out.CodominantCodominant = mySNPLD(geno1,geno2,'codominant-codominant');
                 out.DominantCodominant = mySNPLD(geno1,geno2,'dominant-codominant');
                 out.DominantDominant = mySNPLD(geno1,geno2,'dominant-dominant');
             otherwise
                 return;
         end
end

function out = codetoM(geno)
         n0 = sum(geno==0);
         n2 = sum(geno==2);
         if n0<=n2
             out = geno;
         else
             out = abs(geno-2);
         end
end