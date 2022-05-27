
dataset = 'joinedDatasets';
INPUT_DIR = ['../results/asymmetry/meta_analysis/' dataset '/mean_imputed/not_subsampled/par.mat'];
load(INPUT_DIR)
%%
close all

QQGWAS(PAR.P(:,1))


function QQGWAS(pvals)

pvals = pvals(pvals~=0 & ~isnan(pvals));
n = length(pvals);
obs = -log10(sort(pvals));

% Thin pvals
[~,ia] = unique(round(obs,2)); ia = flip(ia);
ib = find(obs>5); iab = sort(unique([ia;ib])); % No thinning above -log10(P) > 5
obs = obs(iab);
exp = (1:n)'/(n+1);
exp = -log10(exp(iab));

% Plotting expectation line
hold on
mx = exp(1);
x = [0,mx+1];
plot(x,x,'k','LineWidth',5)

% Plotting observations
scatter(exp,obs,25,'r','filled'); hold off
xlabel('Expected -log10(p-values)')
ylabel('Observed -log10(p-values)')
xlim([0 mx+1]);
ylim([0 obs(1)+1]);
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])

end