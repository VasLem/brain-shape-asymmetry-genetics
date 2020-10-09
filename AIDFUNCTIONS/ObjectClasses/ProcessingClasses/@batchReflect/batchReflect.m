classdef batchReflect < batchConvertor
    properties
        Direction = 1;
    end
    methods %Constructor
        function obj = batchReflect(varargin)
          obj = obj@batchConvertor(varargin{:});
          obj.InputFormat = 'Mat Object';
          obj.OutputFormat = 'Mat Object';
          obj.Overwrite = false;
        end
    end
    methods% Interface Functions
        function scan = innerFunction(obj,scan)
            scan.Vertices(obj.Direction,:) = -1*scan.Vertices(obj.Direction,:);
            faces = scan.Faces;
            scan.Faces(2,:) = faces(3,:);
            scan.Faces(3,:) = faces(2,:);
            if ~isempty(scan.PoseLM)      
                scan.PoseLM.Vertices(obj.Direction,:) = -1*scan.PoseLM.Vertices(obj.Direction,:);
                if scan.PoseLM.nrV == 5
                   disp('Mirror 5');
                   scan.PoseLM.Vertices = scan.PoseLM.Vertices(:,[2 1 3 5 4]); 
                elseif scan.PoseLM.nrV == 7
                   disp('Mirror 7'); 
                   scan.PoseLM.Vertices = scan.PoseLM.Vertices(:,[2 1 3 5 4 7 6]);
                elseif scan.PoseLM.nrV == 12
                   disp('Mirror 12'); 
                   scan.PoseLM.Vertices = scan.PoseLM.Vertices(:,[4 3 2 1 7 6 5 10 9 8 11 12]);
                elseif scan.PoseLM.nrV == 14
                   disp('Mirror 14'); 
                   scan.PoseLM.Vertices = scan.PoseLM.Vertices(:,[4 3 2 1 7 6 5 10 9 8 11 12 14 13]);   
                else
                    orig = clone(scan.PoseLM);        
                    [N,D] = knn(kde(scan.PoseLM.Vertices,5),orig.Vertices,1);
                    scan.PoseLM.Vertices = scan.PoseLM.Vertices(:,N);
                end            
            end       
        end     
    end
end