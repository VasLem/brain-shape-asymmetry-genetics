function [downsampledLandmarksIndices, reducedFaces, toUpsampleLandmarksIndices] =  getDownsampledLandmarksIndices(RefScan,reduce, expandSearch)
% Given a Shape3D object (RefScan), *utmost* reduce it with a specific ratio, returning the downsampled indices of the vertices (downsampledLandmarksIndices)
% that can be used to convert the RefScan vertices to the downsampled ones (ie. downsampled = RefScan.Vertices(downsampledLandmarksIndices,:)).
% In addition, a reducedFaces array is returned, that includes the faces that correspond to the reduction applied. Note that the resulting representation will not
% be dense, unless expandSearch is True, which causes a greedy search with a bigger search space, requiring more time and memory, but also reducing the number of
% resulting points. When expandSearch is false, only 1/reduce^2 faces will be scanned as neighbors of each face on the Reference Scan. 
% If it is true, 1/reduce^4 faces are scanned.

disp('..Constructing Search Space..')
bs.faces = RefScan.Faces;
bs.vertices = single(RefScan.Vertices);
disp("...Constructing Centroids..")
tic;
% Apply transformation on the template to receive the undersampled copy as Matlab provides it
rbs = reducepatch(bs.faces,bs.vertices, reduce);
rbs.vertices = single(rbs.vertices);

% Compute centroids of the faces of both bs and rbs
bs.centers = (bs.vertices(bs.faces(:,1),:) + bs.vertices(bs.faces(:,2),:) + bs.vertices(bs.faces(:,3),:))/3;
rbs.centers = (rbs.vertices(rbs.faces(:,1),:) + rbs.vertices(rbs.faces(:,2),:) + rbs.vertices(rbs.faces(:,3),:))/3;
toc;
disp("...Constructing Adjacency..")
tic;
%Compute the  adjacency matrix of the faces. Two faces are adjacent if they share a corner/vertex
adjacency = false(size(bs.faces,1), size(bs.faces,1));
[path,ID] = setupParForProgress(size(rbs.centers, 1));
parfor c=1:size(bs.vertices, 1)
    neighbors = find(sum(bs.faces == c,2) > 0);
    [R, C] = meshgrid(neighbors, neighbors);
    ind = sub2ind(size(adjacency), R, C);
    adjInd{c} = ind;
end
closeParForProgress(path,ID)

for c=1:size(bs.vertices, 1)
   adjacency(adjInd{c}) = true;
end
adjacency = sparse(adjacency);
toc;
disp("...Constructing Path Matrix..")
tic;
% Compute the path matrix for the required reduction ratio. 
% The j-th element of the i-th row is non zero if there is a path of at least 1/reduce steps
% from the i-th face to the j-th face, given the adjacency matrix. 
connectivity = adjacency ^ round(1/reduce) > 0;
% Compute the distance between the centroids of the rbs and the centroids of the bs. Result is a matrix of indices
% where the i-th element of the j-th column is the i-th closest, 
% based on euclidian distance, centroid  in bs to the j-th centroid
if expandSearch
    searchSpace = round(1/reduce^4);
else
    searchSpace = round(1/reduce^2);
end
toc;
disp("...Finding closest centroids..")
tic;
[~,faces_bs2rbs_inds] = pdist2(bs.centers, rbs.centers,"Euclidean",'Smallest',searchSpace); 
toc;
disp('..Finding optimal mapping..')
tic;
% Combine the path matrix with the  last result, keeping only euclidian neighbors of the face X of bs closest to the i-th face of rbs
% that do exist in the path matrix X-th row. After, retrieve the vertices in rbs that are closest to the vertices of the filtered faces.
rbs_vertices_inds_mapped_on_bs = zeros(size(rbs.centers));
[path,ID] = setupParForProgress(size(rbs.centers, 1));
parfor c=1:size(rbs.centers, 1)
    indicesToFilter = faces_bs2rbs_inds(:,c);
    gdists = connectivity(indicesToFilter(1), indicesToFilter(2:end));
    all_faces_corrsepondence_on_bs =  [indicesToFilter(1); indicesToFilter([false,gdists])];
    all_faces_corrsepondence_on_bs_vertices_inds = bs.faces(all_faces_corrsepondence_on_bs,:);
    all_faces_corrsepondence_on_bs_vertices = bs.vertices(all_faces_corrsepondence_on_bs_vertices_inds,:);
    to_match = rbs.vertices(rbs.faces(c,:),:);
    [~, vertices_correspondence_inds] = pdist2(all_faces_corrsepondence_on_bs_vertices, to_match,"Euclidean",'Smallest', 1);
    rbs_vertices_inds_mapped_on_bs(c,:) = all_faces_corrsepondence_on_bs_vertices_inds(vertices_correspondence_inds);
end
closeParForProgress(path,ID)
toc;
disp("..Computing Required Transformation..")
tic;
% Construct transformation indicators
[downsampledLandmarksIndices, ~, toUpsampleLandmarksIndices] = unique(rbs_vertices_inds_mapped_on_bs);
% rbs_vertices_inds_mapped_on_bs_faces 
[reducedFaces, ~] = find(reshape(rbs_vertices_inds_mapped_on_bs,[],1)'==downsampledLandmarksIndices);
reducedFaces = reshape(reducedFaces,size(rbs_vertices_inds_mapped_on_bs));
toc;
end