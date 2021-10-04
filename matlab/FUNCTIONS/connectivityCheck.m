function output = connectivityCheck(Adjacency,LABELS,MASK)
% Adjacency = AdjacencyMatrix;
index = 1:9; 
u_levels = index(end);% 9 levels to analyze
ULABELS = LABELS(index,:);
UHI = HierarchicalInterface;
UHI.nL = u_levels;
UMASK = MASK(1:UHI.nLC);
%disp(['Number of total clusters: ' num2str(length(find(UMASK)))]);
% Check connected component labeling 
AdjacencyMatrix = full(Adjacency); 
output = zeros(length(find(UMASK)),3);
list  = find(UMASK);
   for i=1:length(list)
        [l,c] = Ind2LC(UHI,list(i));
        %disp(['Check Connectivity for Segment for Level ' num2str(l) ' Cluster ' num2str(c)]);
        CLInd = ULABELS(l,:);
        SubInd = find(CLInd==c);
        AdjacencySeg = AdjacencyMatrix(SubInd,SubInd);
        % Connected component labeling: https://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/
        AdjacencySeg(1:1+size(AdjacencySeg,1):end) = 1;  % make diag = 1
        [p,q,r,s] = dmperm(AdjacencySeg);
        numDisjoint = length(r)-1;
        %firstConnectedSet = p(r(1):r(2)-1);
        %secondConnectedSet = p(r(2):r(3)-1);...
        sizeFirstGrouping = r(2)-1;
        percFirstGrouping = sizeFirstGrouping/length(SubInd);
        output(i,1)= list(i);
        output(i,2)= numDisjoint;
        output(i,3)= percFirstGrouping;
   end
end 

