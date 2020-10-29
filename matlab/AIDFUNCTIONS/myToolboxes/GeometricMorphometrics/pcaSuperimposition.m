function Floating = pcaSuperimposition(Target,pcamodel,varargin)
         % read options
         options = readVarargin(varargin{:});
         if isempty(options.index), options.index = 1:Target.nVertices;end
         % create initial result
         Floating = clone(Target);
         if ~(Target.nVertices==pcamodel.Average.nVertices), error('A superimposition cannot be done between shapes that do not share the same amount of vertices'); end       
         % initialize inlier beliefs
         winit = zeros(1,Floating.nVertices);
         winit(options.index)= 1;% activate only points in the index
         w = winit;sumw = sum(w);
         if options.display
            Floating.VertexValue = w;
            v = viewer(Floating); viewer(Target,v);
            Target.ViewMode = 'points';
            Target.SingleColor = [0.8 0.8 0.8];
            Target.ColorMode = 'single';
            Floating.ViewMode = 'Wireframe';
            Floating.ColorMode = 'Indexed';
            pause(1);
         end
         % setting options
         pcamodel.Neta = options.neta;
         pcamodel.Clamp = options.clamp;
         pcamodel.ClampValue = options.clampvalue;
         for iter = 1:1:options.iter% perform ten times
            % get PCA projection
             update = compute_morphable_transformation(pcamodel,Target,w);
             Floating.Vertices = update.Vertices;
            % get remaining distances
             dist = sqrt(sum((Floating.Vertices-Target.Vertices).^2,2))';
            % update sigma of inlier gaussian pdf
             sigma = sqrt(sum(w.*dist.^2)/sumw);
            % update Level outlier uniform pdf  
             L = (1/sqrt(((2*pi)^2)*det(sigma)))*exp(-0.5*options.kappa^2);
            % update inlier beliefs
             ip = normpdf(dist,0,sigma);
            % update outlier beliefs 
             op = repmat(L,1,Floating.nVertices);
            % update weights 
             w = ip./(ip+op);
             w = w.*winit;sumw = sum(w);% incorporate deterministic index by dot multiplication with winit
             if options.display, Floating.VertexValue = w;pause(0.1); end
         end
         Floating.UserData.InlierBeliefs = w;
end

function options = readVarargin(varargin)
         Input = find(strcmpi(varargin, 'options'));
         if ~isempty(Input), options = varargin{Input+1}; return;end    
         Input = find(strcmpi(varargin, 'kappa'));
         if isempty(Input),options.kappa = +inf;else, options.kappa = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'neta'));
         if isempty(Input),options.neta = 1;else, options.neta = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'clamp'));
         if isempty(Input),options.clamp = false;else, options.clamp = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'clampvalue'));
         if isempty(Input),options.clampvalue = 2;else, options.clampvalue = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'display'));
         if isempty(Input),options.display = false;else, options.display = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'index'));
         if isempty(Input),options.index = [];else, options.index = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'iter'));
         if isempty(Input),options.iter = 10;else, options.iter = varargin{Input+1};end
end