function LabelTotal=splitFromSimilarityMatrixv2(SMatrix,n_levels,ADJ)

% This function creates sub-modules of the face through iterative runs of
% the spectral clusteinrg (Ng et al.)
% INPUT
% faces: m x n, m participants, n coordinates of points
% RefScan: 1x1 meshObj

% needs parallel computing on

% INIATIALIZATION

nV = size(SMatrix,1);
LabelTotal = ones(n_levels+1,nV);
n_cluster = 1;
%nV = size(SMatrix,1);
CMatrix = ones(nV,nV);
if nargin<3, ADJ = []; end
connect = true;
if isempty(ADJ), connect = false;end
for level=1:n_levels
    % level = 1;
    % level = 3;
    clustercount = 0;
    UMatrix = SMatrix.*CMatrix;
    for cluster=1:n_cluster
        % cluster = 1;
        a = [datestr(now) ' Level ', num2str(level), ' start', ' cluster ', num2str(cluster), ' start'];
        disp(a)
        pos = find(LabelTotal(level,:)==cluster);
        if LabelTotal(level:end,pos)>0
            if length(pos)>1
               RV = UMatrix(pos,pos);
               label = SpectralAffinity(double(RV)); 
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
    n_cluster = clustercount;
    % UPDATING CONNECTIVITY MATRIX
    if connect&&(level<n_levels)
        % NEED TO UPDATE CONNECTIVITY MATRIX
        disp('Updating Connectivity');
        values = LabelTotal(level+1,:);
        ModMatrix = eye(nV,nV);
        for i=1:nV
           ModMatrix(i,:) = values==values(i); 
        end
        clear values;
        D = pairWiseGeodesicDistancesFromAdjacency(ADJ.*ModMatrix);
        %D = FastFloyd(1./RefScan.Adjacency);
        CMatrix = single(D>=0);
        CMatrix(CMatrix==0) = 0.01;
        %CMatrix(CMatrix>1) = 1;
        clear D;
    end
    %waitbar(level/n_levels,h);
end
%close(h)



