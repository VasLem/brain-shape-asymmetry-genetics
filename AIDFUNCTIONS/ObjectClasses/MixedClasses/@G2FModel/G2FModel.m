classdef G2FModel < superClass
   properties
      Model;
      Shape;
      Sex;
      RIPS;
      Anc;
      RIPA;
      GNAMES;
      G;
      RIPG;
      BF = 1.5;
   end
   properties (Dependent = true)
      Findex;
      Mindex;
      AncP;
      AncS;
      nrG;
      nrS;
   end
   methods % Constructor
        function obj = G2FModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % special setting and getting
       function out = get.Findex(obj)
           out = find(obj.Sex==1);
       end
       function out = get.Mindex(obj)
           out = find(obj.Sex==-1);
       end
       function out = get.AncP(obj)
          [out,~] = polyfit(obj.Anc,obj.RIPA,1);
       end
       function out = get.AncS(obj)
          [~,out] = polyfit(obj.Anc,obj.RIPA,1);
       end
       function out = get.nrG(obj)
           out = length(obj.GNAMES);
       end
       function out = get.nrS(obj)
          out = length(obj.Sex); 
       end
   end
   methods % Interface function
       
   end
end