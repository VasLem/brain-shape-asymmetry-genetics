function G2FMontage(obj,pred)
disp('-------------------JANUS FACIAL DNA Montage-------------------');
disp('');
disp('Created by Peter Claes and Mark Shriver, Copyright 2013');
disp('--------------------------------------------------------------');
%% retrieve facial sex
if strcmp(pred.Sex,'M')
   index = obj.Mindex;
elseif strcmp(pred.Sex,'F')
   index = obj.Findex;
else
   error('Gender not recognized, please enter M or F only');
end
pred.RIPS = mean(obj.RIPS(index));
pred.RIPSE = std(obj.RIPS(index));
disp(['Given Sex  = ' pred.Sex]);
disp(['Estimated Facial Sex = ' num2str(pred.RIPS)]);
disp(['Estimated error on Facial Sex = ' num2str(pred.RIPSE)]);
disp('--------------------------------------------------------------'); 
%% retrieve facial ancestry         
[pred.RIPA,pred.RIPAE] = polyval(obj.AncP,pred.Anc,obj.AncS);
disp(['Given Genomic Ancestry = ' num2str(pred.Anc)]);
disp(['Estimated Facial Ancestry = ' num2str(pred.RIPA)]);
disp(['Estimated error on Facial Ancestry = ' num2str(pred.RIPAE)]);
disp('--------------------------------------------------------------'); 
%% Make the base face
getBaseFace(obj,pred);
disp('Base Face generated...');
%% retrieve facial gene effects
disp('Processing GENO TYPES');
f = statusbar('Processing...');drawnow;
GeneEffects = cell(1,pred.nrG);
good = [];
for g=1:pred.nrG
    if pred.G(g)==0, continue; end
    i = find(strcmp(pred.GNAMES{g},obj.GNAMES));
    if isempty(i), continue; end
    G.Name = pred.GNAMES{g};
    G.Type = pred.G(g);
    A = [obj.RIPS obj.RIPA obj.RIPG(:,i)];
    B = obj.Shape;
    [A,B] = eliminateNAN(A,B);
    G.Effect = createGeneEffects(A,B,pred.BaseFace,obj.BF,[pred.RIPS pred.RIPA]);
    switch G.Type
        case -1
            G.rightobj = G.Effect.Morph1;
            G.rightobj.Tag = [G.Name '_1'];         
            G.wrongobj = G.Effect.Morph2;
            G.wrongobj.Tag = [G.Name '_2'];
        case 1
            G.rightobj = G.Effect.Morph2;
            G.rightobj.Tag = [G.Name '_2'];
            G.wrongobj = G.Effect.Morph1;
            G.wrongobj.Tag = [G.Name '_1'];
    end
    GeneEffects{g} = G;
    good = [good g]; %#ok<AGROW>
    statusbar(g/pred.nrG,f);drawnow;
end
delete(f);
%% filtering missed gene effects
GeneEffects = GeneEffects(good);
nrG = length(GeneEffects);
disp('--------------------------------------------------------------');
%% Gene effect overlaying
disp('Overlaying Gene Effects...');
shapeR = zeros(pred.BaseFace.nrV*3,nrG);
weightsR = zeros(pred.BaseFace.nrV*3,nrG);
for i=1:1:nrG
    shapeR(:,i) = GeneEffects{i}.rightobj.Vertices(:);
    tmp = repmat(GeneEffects{i}.Effect.STATS.LocalS,3,1);
    weightsR(:,i) = tmp(:);
end
pred.PredFace = clone(pred.BaseFace);
shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
pred.PredFace.Vertices = reshape(shaperes,3,pred.BaseFace.nrV);
pred.PredFace.Value = vDistances(pred.PredFace,pred.BaseFace);
% disp('Amplifying Identity...');
% out = amplifyIdentity(out,AF);
% if ass
%    disp('Assessing Identity...');
%    out = assessIdentity(out);
% end
disp('--------------------------------------------------------------');
disp('DONE');
end

function out = createGeneEffects(A,shape,RefScan,BF,CondValues)
         % multiple regression to obtain multiple R-Squared
            [~,~,~,~,~,pctvar] = plsregress(A,shape,size(A,2));
            R2Tot = sum(pctvar(2,:));
            STATS.R2Tot = R2Tot;
         % Conditioning regression to obtain reference face according to CondValues
            AC = A(:,1:end-1);
            [~,~,~,~,MC,pctvarC,~,statsC] = plsregress(AC,shape,size(AC,2));
            E = statsC.Yresiduals;
            mC = mean(AC);
            mShape = mean(shape);
            if isempty(CondValues)
               CondValues = mC; 
            end
            DX = mC-CondValues;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            scanC = clone(RefScan);
            delete(scanC.PoseLM);
            scanC.Vertices = reshape(NY,3,length(NY)/3);
            R2Red = sum(pctvarC(2,:));
            R2Add = R2Tot-R2Red;
            STATS.R2Add = R2Add;
         % reduced model construction   
            AR = A(:,end);
            [~,~,~,~,~,~,~,statsR] = plsregress(AC,AR,1);
            AR = statsR.Yresiduals;
            [~,~,~,~,MR,pctvar,~,stats] = plsregress(AR,E,1);
            R2Part = pctvar(2);
            STATS.R2Part = R2Part;
            [LocalE,LocalS] = getLocalStatistic(MR,stats.Yresiduals,E,'shape');
            STATS.LocalE = LocalE;
            STATS.LocalS = LocalS;
          % creating morphs in reduced model space   
            Y = mean(E);
            m = mean(AR);mv = std(AR);
            range(1) = BF*(m-3*mv);
            range(2) = BF*(m+3*mv);
            counter = 0;
            for k=1:1:2
                V = range(k);
                counter = counter+1;
                DX = m-V;
                DY = DX*MR(2:end,:);
                NY = Y-DY;
                scan = clone(RefScan);
                delete(scan.PoseLM);
                scan.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
                eval(['morph' num2str(counter) ' = clone(scan);']);
            end
            out.Morph1 = morph1;
            out.Morph2 = morph2;
            out.STATS = STATS;
end

function [LocalE,LocalS] = getLocalStatistic(M,R,B,type)
         [n,nB] = size(B); 
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         switch lower(type)
             case 'value'
                 LocalE = M(2,:);
                 LocalS = SSR./SST;
              case 'shape'
                 LocalE = zeros(4,nB/3);
                 LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
                 LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
                 SSR = reshape(SSR,3,nB/3);
                 SSR = sum(SSR);
                 SST = reshape(SST,3,nB/3);
                 SST = sum(SST);
                 LocalS = SSR./SST;
             otherwise
                 error('unknown type');
         end        

end