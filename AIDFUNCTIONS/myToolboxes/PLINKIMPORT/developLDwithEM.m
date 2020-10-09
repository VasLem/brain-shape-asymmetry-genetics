% develop LD with EM
load('D:\Dropbox\mySHARE\NEW\AUGUSTUPDATE\PHENOTYPES\ExampleGenotypes.mat');

out = mySNPLD(geno1,geno2);


new1 = geno1;
new2 = geno1;

index = randperm(length(new1));
out = mySNPLD(new1,new2(index));
