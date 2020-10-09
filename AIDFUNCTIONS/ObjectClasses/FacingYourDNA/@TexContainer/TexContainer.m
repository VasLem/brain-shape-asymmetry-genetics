classdef TexContainer < superClass
   % A container class containing all the parameter values
   properties
      TrackID = 'X';
      TexSpace = [];
   end
   properties (Dependent = true, Hidden = true)
      DepVar;
      RefScan;
      AvgTex;
      UV;
      n;
   end
   properties (Hidden = true)
       
   end
   methods % Constructor
        function obj = TexContainer(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % GENERAL GETTING\SETTING
       function out = get.TexSpace(obj)
           out = obj.TexSpace;
           if ~superClass.isH(out), out = []; end
       end       
       function out = get.n(obj)
          if isempty(obj.TexSpace), out = 0; return; end
          out = obj.TexSpace.n;
       end
       function out = get.DepVar(obj)
            if isempty(obj.TexSpace), out = [];return; end
            out = obj.TexSpace.Tcoeff./repmat(obj.TexSpace.EigStd',obj.n,1);
       end
       function out = get.RefScan(obj)
            if isempty(obj.TexSpace), out = []; return; end
            out = obj.TexSpace.RefScan;
       end
       function out = get.UV(obj)
           if isempty(obj.RefScan), out = []; return; end
           out = obj.RefScan.UV;
       end
       function out = get.AvgTex(obj)
           if isempty(obj.TexSpace), out = []; return; end
           out = obj.TexSpace.AvgTexture;
       end
   end
   methods % INTERFACING FUNCTIONS
       function out = extractTexture(obj,A,X)
                avgA = nanmean(A);B = obj.DepVar;
                X(isnan(X)) = avgA(isnan(X));
                [A,B] = TexContainer.eliminateNAN(A,B);
                [~,~,~,~,M] = plsregress(A,B,size(A,2));
                C = [1 X]*M;
                out = getScan(obj.TexSpace,C.*obj.TexSpace.EigStd');
                out.ColorMode = 'Texture';
       end
   end
   methods (Static = true)
      function [A,B] = eliminateNAN(A,B)
         index = (1:size(A,1));
         [i,~] = find(isnan(A));
         i = unique(i);
         index = setdiff(index,i);
         A = A(index,:);
         B = B(index,:);
      end 
   end
end