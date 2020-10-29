classdef PredModelContainer < superClass
   % A container class containing all the parameter values
   properties
       TrackID = 'X';
   end
   properties (Hidden = true)
      BaseCont = [];
      SNPCont = [];
      SNPBaseCont =  [];
      TexCont = [];
      SameBase = 0;
   end
   methods % Constructor
        function obj = PredModelContainer(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % GETTING/SETTING
       function out = get.BaseCont(obj)
           out = obj.BaseCont;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.SNPCont(obj)
           out = obj.SNPCont;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.TexCont(obj)
           out = obj.TexCont;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.SNPBaseCont(obj)
           if obj.SameBase, out = obj.BaseCont; return; end
           out = obj.SNPBaseCont;
           if ~superClass.isH(out), out = []; end
       end
   end
end