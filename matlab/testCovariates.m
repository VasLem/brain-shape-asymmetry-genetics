DATA_DIR = '../SAMPLE_DATA/';
UKBIOBANK = 'UKBIOBANK';
COV_GENO_PATH = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
cov_stage00data = load(COV_GENO_PATH);

UKBIOBANK = 'MY_UKBIOBANK';
COV_GENO_PATH = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];
cov_batch2_2021_data = load(COV_GENO_PATH);
%%
table1 = array2table(cov_stage00data.COV.DATA,"VariableNames",cellfun(@(x)x(1:min(length(x),40)),cov_stage00data.COV.Names, 'UniformOutput',false));
%%
table2 = array2table(cov_batch2_2021_data.COV.DATA,"VariableNames",cellfun(@(x)x(1:min(length(x),40)),cov_batch2_2021_data.COV.Names, 'UniformOutput',false));
%%
figure;
subplot(2,1,1);
hist(table1.Height)
subplot(2,1,2);
hist(table2.("Standing height"))
%%
close all
drawhist(table1.("KUL GENO PC 1"), table2.("KUL GENO PC 1"))
%%
drawhist(table1.("Genetic sex"), table2.("Genetic sex"))
%%
function drawhist(H1, H2)
map = brewermap(2,'Set1'); 
figure
histogram(H1,min(min(H1),min(H2)):.01:max(max(H1),max(H2)),'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
hold on
histogram(H2,min(min(H1),min(H2)):.01:max(max(H1),max(H2)),'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
box off
axis tight
legend("1","2",'location','northwest')
legend boxoff
end
