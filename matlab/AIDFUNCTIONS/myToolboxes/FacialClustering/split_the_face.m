function LabelTotal=split_the_face(faces,avScan,n_levels,superimpose)

% This function creates sub-modules of the face through iterative runs of
% the spectral clusteinrg (Ng et al.)
% INPUT
% faces: m x n, m participants, n coordinates of points
% RefScan: 1x1 meshObj

% needs parallel computing on
avFace = avScan.Vertices;
[nObs,nV] = size(faces);
nbLM = size(avFace,1);

Mfaces = reshape(faces',3,nbLM,nObs);

%h=waitbar(0,'please wait');
% INIATIALIZATION
LabelTotal = ones(n_levels+1,7150);
n_cluster = 1;

%RefScanVerticesReshape = avFace(:);
%Displacement = bsxfun(@minus, repmat(RefScanVerticesReshape',nObs,1), faces);
%RVMatrix= rand(nbLM,nbLM);
RVMatrix=buildRVmatrix(faces);
CMatrix = ones(nbLM,nbLM);
for level=1:n_levels
    % level = 1;
    % level = 2;
    clustercount = 0;
    RVMatrix = RVMatrix.*CMatrix;
    for cluster=1:n_cluster
        % cluster = 1;
        a = [datestr(now) ' Level ', num2str(level), ' start', ' cluster ', num2str(cluster), ' start'];
        disp(a)
        pos = find(LabelTotal(level,:)==cluster);
        if LabelTotal(level:end,pos)>0
            if length(pos)>1
                if superimpose
                   %temp = repmat((pos-1)*3+1,3,1)+repmat([0;1;2],1,length(pos));
                   if level==1
                      RV = RVMatrix;
                   else
                      scan = crop(avScan,'VertexIndex',pos); 
                      Mredfaces = Mfaces(:,pos,:);
                      Mredfaces = reshape(Mredfaces,length(pos)*3,nObs);
                      SymSpace = shapePCA;
                      SymSpace.RefScan = clone(scan);
                      getAverage(SymSpace,Mredfaces);
                      D = LSGenProcrustes(SymSpace,Mredfaces,true,3,scan);
                      RV=buildRVmatrix(D');
                   end
                else 
                  RV = RVMatrix(pos,pos);
                end
                
                label = SpectralAffinity(RV); 
                label(label==1)=clustercount+1;
                label(label==2)=clustercount+2;
                clustercount = clustercount + 2;
                LabelTotal(level+1,pos)=label;
            else
                LabelTotal(level+1:end,pos)=-cluster/(level);
            end
        end
        a = [datestr(now) ' Level ', num2str(level), ' done', ' cluster ', num2str(cluster), ' done'];
        disp(a)
    end
    % UPDATING CONNECTIVITY MATRIX
    
    n_cluster = clustercount;
    %waitbar(level/n_levels,h);
end
%close(h)



