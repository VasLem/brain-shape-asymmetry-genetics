classdef scaledRigidTM < TM
    properties %(GetAccess = private, SetAccess = private)
        Rotation = eye(3,3);
        Translation = zeros(3,1);
        Scale = 1;
        d = [0.001; 0.001; 0.001; 0.1; 0.1; 0.1;0.01];% Parameter deviations to compute numerical gradients       
        c = [];
        Fields = {'Vertices'};% rigidTM only affects the Vertices of LMObj and meshObj
    end
    properties (Dependent = true)
        P;% Transformation parameters (three rotation angles and three translations
        Rx;
        Ry;
        Rz;
        Tx;
        Ty;
        Tz;
        HomMatrix;
    end
    methods %Constructor
        function obj = scaledRigidTM(varargin)
          obj = obj@TM(varargin{:});
          if nargin > 0
            Input = find(strcmp(varargin, 'c'));if ~isempty(Input), obj.c = varargin{Input+1}; end
            Input = find(strcmp(varargin, 'Rotation'));if ~isempty(Input), obj.Rotation = varargin{Input+1}; end
            Input = find(strcmp(varargin, 'Translation'));if ~isempty(Input), obj.Translation = varargin{Input+1}; end
            Input = find(strcmp(varargin, 'Scale'));if ~isempty(Input), obj.Scale = varargin{Input+1}; end
          end
        end
    end
    methods % Setting and Getting
        function out = get.P(obj)
                 out = zeros(7,1);
                 out(1:3) = rodrigues(obj.Rotation);
                 out(4:6) = obj.Translation;
                 out(7) = obj.Scale;
        end
        function obj = set.P(obj,in)
                 if size(in,1)==1, in = in'; end
                 obj.Rotation = rodrigues(in(1:3));
                 obj.Translation = in(4:6);
                 obj.Scale = in(7);
        end
        function r = get.Rx(obj)
            r = obj.P(1);
        end
        function obj = set.Rx(obj,r)
            obj.P(1) = r;
        end
        function r = get.Ry(obj)
            r = obj.P(2);
        end
        function obj = set.Ry(obj,r)
            obj.P(2) = r;
        end
        function r = get.Rz(obj)
            r = obj.P(3);
        end
        function obj = set.Rz(obj,r)
            obj.P(3) = r;
        end
        function t = get.Tx(obj)
            t = obj.P(4);
        end
        function obj = set.Tx(obj,t)
            obj.P(4) = t;
        end
        function t = get.Ty(obj)
            t = obj.P(5);
        end
        function obj = set.Ty(obj,t)
            obj.P(5) = t;
        end
        function t = get.Tz(obj)
            t = obj.P(6);
        end
        function obj = set.Tz(obj,t)
            obj.P(6) = t;
        end
        function out = get.HomMatrix(obj)
            out = eye(4,4);
            out(1:3,1:3) = obj.Scale*obj.Rotation;
            out(1:3,4) = obj.Scale*obj.Translation;
        end
        function obj = set.HomMatrix(obj,in)
            obj.Scale = 1;
            obj.Rotation = in(1:3,1:3);
            obj.Translation = in(1:3,4);
        end
    end 
    methods % InterFace functions
        function copy(obj,cobj)
                 copy@TM(obj,cobj);
                 cobj.c = obj.c;
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@TM(obj);
            struc.c = obj.c;
        end
        function obj = struc2obj(obj,struc)
                 obj = struc2obj@TM(obj,struc);
                 obj.c = struc.c;
        end
    end
end % classdef

%         function obj = set.c(obj,p)
%                  p = TM.getPoints(p);
%                  obj.c = mean(p,2);
%                  %obj.c = [0;0;0];
%         end