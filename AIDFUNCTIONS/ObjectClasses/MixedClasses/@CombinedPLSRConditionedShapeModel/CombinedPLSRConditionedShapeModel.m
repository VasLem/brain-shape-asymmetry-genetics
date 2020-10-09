classdef CombinedPLSRConditionedShapeModel < PLSRSuper
   properties (Dependent = true)
      Model; 
   end
   properties
       CYModel = PLSRShapeModel;
       XModels = {};
       ActiveX= 0;
   end
   properties (Dependent = true)
      nrXModels;
      X;
      XNames;
      Y;
      C;
      CNames;
      CAvg;
      CStd;
      nrC;
      Average;
   end
   methods % Constructor
        function obj = CombinedPLSRConditionedShapeModel(varargin)
            obj = obj@PLSRSuper(varargin{:});
        end
   end
   methods % special setting and getting
        function out = get.CYModel(obj)
           out = obj.CYModel;
           if ~superClass.isH(out)
               new = PLSRShapeModel;
               obj.CYModel = new;
               out = new;
               if obj.nrXModels == 0, return; end
               for i=1:1:obj.nrXModels
                   obj.XModels{i}.CYModel = new;
               end
           end
         end
        function out = get.Model(obj)
           if isempty(obj.CYModel), out = []; return; end
           out = obj.CYModel.Model;
           if ~superClass.isH(out), out = []; end
        end
        function obj = set.Model(obj,in)
                 obj.CYModel.Model = in;
        end
        function out = get.nrXModels(obj)
           out = length(obj.XModels); 
        end
        function obj = set.X(obj,in)
           nrIn = size(in,2);
           obj.XModels = cell(1,nrIn);
           for i=1:1:nrIn
%                obj.XModels{i} = PLSRConditionedShapeModel;
%                obj.XModels{i}.CYModel = obj.CYModel;
%                obj.XModels{i}.X = in(:,i);
%                obj.XModels{i}.CXModel.Parent = obj.XModels{i};
               new = PLSRConditionedShapeModel;
               new.CYModel = obj.CYModel;
               new.X = in(:,i);
               %new.CXModel.Parent = new;
               obj.XModels{i} = new;
           end
        end
        function out = get.X(obj)
           if obj.nrXModels==0, out = []; return; end
           out = zeros(obj.nrO,obj.nrXModels);
           for i=1:1:obj.nrXModels
              out(:,i) = obj.XModels{i}.X;  
           end
        end
        function obj = set.XNames(obj,in)
           if obj.nrXModels==0
              nrIn = size(in,2);
              obj.XModels = cell(1,nrIn);
               for i=1:1:nrIn
%                    obj.XModels{i} = PLSRConditionedShapeModel;
%                    obj.XModels{i}.CYModel = obj.CYModel;
%                    obj.XModels{i}.XNames = in(i);
%                    obj.XModels{i}.CXModel.Parent = obj.XModels{i};
                   new = PLSRConditionedShapeModel;
                   new.CYModel = obj.CYModel;
                   new.XNames = in(i);
                   %new.CXModel.Parent = new;
                   obj.XModels{i} = new;
               end
               return;
           end
           for i=1:1:obj.nrXModels
               obj.XModels{i}.XNames = in(i);
           end
        end
        function out = get.XNames(obj)
           if obj.nrXModels == 0, out = []; return; end
           out = cell(1,obj.nrXModels);
           for i=1:1:obj.nrXModels
               tmp = obj.XModels{i}.XNames;
               if isempty(tmp), continue; end
               out{i} = obj.XModels{i}.XNames{:}; 
           end
        end
        function out = get.Y(obj)
            out = obj.CYModel.Y;
        end
        function out = get.Average(obj)
           out = obj.CYModel.Average;
        end
        function obj = set.C(obj,in)
            obj.CYModel.X = in;
        end
        function out = get.C(obj)
            out = obj.CYModel.X;
        end
        function obj = set.CNames(obj,in)
            obj.CYModel.XNames = in;
        end
        function out = get.CNames(obj)
            out = obj.CYModel.XNames;
        end
        function out = get.CAvg(obj)
           out = obj.CYModel.XAvg; 
        end
        function out = get.CStd(obj)
           out = obj.CYModel.XStd;
        end
        function out = get.nrC(obj)
            out = obj.CYModel.nrX;
        end
   end
   methods % Interface functions
       function update(obj)
           updateCYModel(obj);
           updateXModels(obj);
       end
       function updateCYModel(obj)           
            update(obj.CYModel);
       end
       function updateXModels(obj)
          if obj.nrXModels==0, return; end
          for i=1:1:obj.nrXModels
             updateX(obj.XModels{i}); 
          end
       end
       function [out,scan] = dC(obj,in,Y,index)
                if nargin<4, index = (1:obj.nrC);end
                [out,scan] = dX(obj.CYModel,in,Y,index); 
       end
       function [out,scan] = changeCY(obj,in,C,Y)
                [out,scan] = changeX(obj.CYModel,in,C,Y);
       end
       function [out] = changeCX(obj,in,C,X)
           if obj.nrXModels == 0, out = []; return; end
           out = zeros(1,obj.nrXModels);
           for i=1:1:obj.nrXModels
            out(i) = changeX(obj.CXModel{i},in,C,X(i));
           end
       end
       function [outY,outX,scan] = changeC(obj,in,C,Y,X)
                [outY,scan] = changeCY(obj,in,C,Y);
                [outX] = changeCX(obj,in,C,X);
       end
       function out = CfromY(obj,in)
                out = XfromY(obj.CYModel,in);
       end
       function [out,scan] = YfromC(obj,in)
                [out,scan] = YfromX(obj.CYModel,in);
       end
       function out = XfromC(obj,in)
           if obj.nrXModels==0, out = []; return; end
           out = zeros(1,obj.nrXModels);
           for i=1:1:obj.nrXModels
               out(i) = XfromC(obj.XModels{i},in);
           end
       end
       function out = CfromX(obj,in) %#ok<*STOUT,*INUSD,*MANU>
           %out = XfromY(obj.CXModel,in); 
           error('This action cannot be performed at this stage');
       end
       function out = XfromYC(obj,Y,C)
           if obj.nrXModels==0, out = []; return; end
           out = zeros(1,obj.nrXModels);
           for i=1:1:obj.nrXModels
               out(i) = XfromYC(obj.XModels{i},Y,C);
           end
       end
       function [out,scan] = YfromXC(obj,X,C)
            error('This action cannot be performed at this stage');
       end
       function out = YResfromC(obj,in,C)
           out = YResfromX(obj.CYModel,in,C);
       end
       function out = XResfromC(obj,in,C)
           if obj.nrXModels==0, out = []; return; end
           out = zeros(1,obj.nrXModels);
           for i=1:1:obj.nrXModels
               out(i) = XResfromC(obj.XModels{i},in,C);
           end
       end
       function [out,scan] = changeX(obj,in,X,Y,C)
           if obj.nrXModels == 0, out = []; scan = []; return; end
            % get Y residual from C
            ybase = YfromC(obj,C);
            ycurres = Y-ybase;
            yres = zeros(obj.nrY,obj.nrXModels);
            for i=1:1:obj.nrXModels
                % get X residual current and in from C
                    xres = XResfromC(obj.XModels{i},in(i),C);
                    xcurres = XResfromC(obj.XModels{i},X(i),C);
                % change X current to X in in reduced model space and create out  
                    yres(:,i) = changeX(obj.XModels{i}.XYModel,xres,xcurres,ycurres);
            end
            yres = sum(yres,2);
            out = ybase + yres';
            if nargout==1, return; end
            scan = getScan(obj.Model,out);
       end
       function [out,scan] = changeSingleX(obj,in,x,Y,C)
           if obj.nrXModels == 0, out = []; scan = []; return; end
           % get Y residual from C
           ybase = YfromC(obj,C);
           ycurres = Y-ybase;
           % get X residual current and in from C
           xres = XResfromC(obj.XModels{obj.ActiveX},in,C);
           xcurres = XResfromC(obj.XModels{obj.ActiveX},x,C);
           % change X current to X in in reduced model space and create out  
           yres = changeX(obj.XModels{obj.ActiveX}.XYModel,xres,xcurres,ycurres);
           out = ybase + yres;
           if nargout==1, return; end
           scan = getScan(obj.Model,out);   
       end
       function reconnectChildren(obj)
           if obj.nrXModels == 0, return; end
           for i=1:1:obj.nrXModels
              reconnectChildren(obj.XModels{i}); 
           end
       end
   end
        
    
    
    
end