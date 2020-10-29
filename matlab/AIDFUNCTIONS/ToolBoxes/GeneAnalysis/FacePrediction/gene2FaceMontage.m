function out = gene2FaceMontage(input,G2FModel,BF,AF,ass)
disp('-------------------JANUS FACIAL DNA Montage-------------------');
disp('');
disp('Created by Peter Claes and Mark Shriver, Copyright 2013');
disp('--------------------------------------------------------------');
out.input = input;
%% retrieve facial sex
if strcmp(input.Sex,'M')
   index = G2FModel.Mindex;
elseif strcmp(input.Sex,'F')
   index = G2FModel.Findex;
else
   error('Gender not recognized, please enter M or F only');
end
out.Fsex = mean(G2FModel.RIPS(index));
out.stdFsex = std(G2FModel.RIPS(index));
disp(['Given Sex  = ' input.Sex]);
disp(['Estimated Facial Sex = ' num2str(out.Fsex)]);
disp(['Estimated error on Facial Sex = ' num2str(out.stdFsex)]);
disp('--------------------------------------------------------------'); 
%% retrieve facial ancestry         
[out.Fanc,out.stdFanc] = polyval(G2FModel.AncP,input.Anc,G2FModel.AncS);
disp(['Given Genomic Ancestry = ' num2str(input.Anc)]);
disp(['Estimated Facial Ancestry = ' num2str(out.Fanc)]);
disp(['Estimated error on Facial Ancestry = ' num2str(out.stdFanc)]);
disp('--------------------------------------------------------------'); 
%% Make the base face
out.BaseFace = getBaseFace([G2FModel.RIPS G2FModel.RIPA],[out.Fsex out.Fanc],G2FModel.Model.Tcoeff,G2FModel.Model);
out.BaseFace.SingleColor = [0.8 0.8 0.8];
%% retrieve facial gene effects
disp('Processing GENO TYPES');
f = statusbar('Processing...');drawnow;
out.GeneEffects = cell(1,input.nrGenes); 
for g=1:input.nrGenes
    G.Gname = input.Genes{g};
    G.Type = input.GenoTypes(g);
    i = find(strcmp(input.Genes{g},G2FModel.Genes));
    A = [G2FModel.RIPS G2FModel.RIPA G2FModel.RIPG(:,i)];
    G.Effect = createGeneEffects(A,G2FModel.Shape,out.BaseFace,BF,[out.Fsex out.Fanc]);
    if G.Type==0, G.Type = G2FModel.HETERO(i); end
    switch G.Type
        case -1
            G.rightobj = G.Effect.Morph1;
            G.rightobj.Tag = [G.Gname '_1'];         
            G.wrongobj = G.Effect.Morph2;
            G.wrongobj.Tag = [G.Gname '_2'];
        case 1
            G.rightobj = G.Effect.Morph2;
            G.rightobj.Tag = [G.Gname '_2'];
            G.wrongobj = G.Effect.Morph1;
            G.wrongobj.Tag = [G.Gname '_1'];
    end
    out.GeneEffects{g} = G;
    statusbar(g/input.nrGenes,f);drawnow;
end
delete(f);
disp('--------------------------------------------------------------');
%% Gene effect overlaying
disp('Overlaying Gene Effects...');
shapeR = zeros(out.BaseFace.nrV*3,input.nrGenes);
weightsR = zeros(out.BaseFace.nrV*3,input.nrGenes);
for i=1:1:input.nrGenes
    shapeR(:,i) = out.GeneEffects{i}.rightobj.Vertices(:);
    tmp = repmat(out.GeneEffects{i}.Effect.STATS.LocalS,3,1);
    weightsR(:,i) = tmp(:);
end
out.PredFace = clone(out.BaseFace);
shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
out.PredFace.Vertices = reshape(shaperes,3,out.BaseFace.nrV);
out.PredFace.Value = vDistances(out.PredFace,out.BaseFace);
disp('Amplifying Identity...');
out = amplifyIdentity(out,AF);
if ass
   disp('Assessing Identity...');
   out = assessIdentity(out);
end
disp('--------------------------------------------------------------');
disp('DONE');
end