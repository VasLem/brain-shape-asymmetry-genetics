function out = myEnsembleLD(rs1,rs2,pop)
         if nargin<3, pop = 'CEU'; end% default population is Great Brittian as reference to European GWAS
         server = 'https://rest.ensembl.org';
         % check if rs or position is given
         if contains(rs1,'rs')
           ext =  ['/ld/human/pairwise/' rs1 '/' rs2 '?population_name=1000GENOMES:phase_3:' pop];
         else
           ext = ['/ld/human/region/' rs1 '..' rs2(3:end) '/1000GENOMES:phase_3:' pop '?'];
         end
         out = webread([server ext]);
end


% rs1 = 'rs6792369';
% rs2 = 'rs1042779';