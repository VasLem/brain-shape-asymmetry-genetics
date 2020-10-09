function [shapeSpace,meshSpace] = importOxfordModel(filename)

    [~,~,raw] = xlsread(filename);

    counter = 3;
    nPC = raw{counter,1};
    counter = counter+1;
    EigVal = cell2mat(raw(counter:counter+nPC-1,1));
    counter = counter+nPC;
    nV = raw{counter,1}; 
    counter = counter+2;
    EigVec = cell2mat(raw(counter:counter+nV*nPC-1,:));
    counter = counter+nV*nPC;
    nT = raw{counter,1};
    counter = counter+1+nT;
    nFaces = raw{counter,1};
    counter = counter+1;
    Faces = cell2mat(raw(counter:counter+nFaces-1,:))+1;
    counter = counter+nFaces;
    nVertices = raw{counter,1};
    counter = counter+1;
    Vertices = cell2mat(raw(counter:counter+nVertices-1,:));

    oxTemplate = shape3D;
    oxTemplate.Vertices = Vertices;
    oxTemplate.Faces = Faces;

    EigVec = EigVec';
    EigVec = EigVec(:);
    EigVec = reshape(EigVec,nV*3,nPC);

    % creating old software PCA object
    RefScan = meshObj;
    RefScan.Vertices = oxTemplate.Vertices';
    RefScan.Faces = oxTemplate.Faces';

    meshSpace = shapePCA;
    meshSpace.RefScan = clone(RefScan);

    meshSpace.EigVal = EigVal;
    meshSpace.EigVec = EigVec;
    
    % creating new software PCA object
    shapeSpace = morphableShape3D;
    shapeSpace.Average = clone(oxTemplate);
    shapeSpace.EigVec = EigVec;
    shapeSpace.EigVal = EigVal;
end