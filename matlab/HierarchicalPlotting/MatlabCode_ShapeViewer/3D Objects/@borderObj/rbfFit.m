function out = rbfFit(obj,varargin)
    % border is the border structure of a facial mesh created with
    % detect_border and DistanceToBorder
    % nrVertices is the amount of non border mesh points you want to use to create the border_rbf
    % default = 3000
    % The total amount of points is increased with the border points
    [nr,acc] = readVarargin(varargin{:});
    if nr<obj.Parent.nrV
       factor = floor(obj.Parent.nrV/nr);
    else
       factor = 1;
    end
    distances = obj.Distance2Border;
    index = (1:factor:obj.Parent.nrV);
    index = union(index,find(distances==0));
    maxdist = max(distances);
    distances=log(1+100*distances/maxdist)/log(101);
    dens.Location = obj.Parent.Location(:,index);
    dens.Value = distances(index);
    %figure;fastrbf_view(dens);
    %dens = uniquePL(dens);
    %rbf = fastrbf_fit(dens,acc,'reduce','messages',0);
    tmprbf = fastRBF;
    tmprbf.FitAccuracy = acc;
    tmprbf.Reduce = true;
    tmprbf.Rho = 0;
    fit(tmprbf,dens);
    if nargout == 1
       out{1} = tmprbf;
    else
       obj.RBF = tmprbf;
    end
    
end

function [nr,acc] = readVarargin(varargin)
    Input = find(strcmp(varargin,'nrVertices'));
    if ~isempty(Input)
       nr = varargin{Input+1};
    else 
       nr = 500;
    end
    Input = find(strcmp(varargin,'Accuracy'));
    if ~isempty(Input)
      acc = varargin{Input+1};
    else 
      acc = 0.1;
    end
end