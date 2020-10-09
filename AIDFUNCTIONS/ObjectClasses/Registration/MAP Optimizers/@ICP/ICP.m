classdef ICP < mapOPT
    % This is the abstract interface class for MAP function optimizers
    % based on a Iterative Corresponding Point principle during Mstep of
    % the optimization
   methods %Constructor
        function obj = ICP(varargin)
          obj = obj@mapOPT(varargin{:});
        end
   end
   methods % InterFace functions
       function out = initialize(obj)
           if nargout == 1
                    obj = clone(obj);
                    out = obj;
           end
           initialize@mapOPT(obj);
           obj.ObjFun.Correspondence = obj.ObjFun.Floating;
           if obj.Show% Showing Correspondences
              obj.ObjFun.Correspondence.ViewMode = 'Points';
              obj.ObjFun.Correspondence.SingleColor = [1 1 1];
              obj.ObjFun.Correspondence.ColorMode = 'Single';
              obj.ObjFun.Correspondence.Axes = obj.ObjFun.Solution.Axes;
              obj.ObjFun.Correspondence.Visible = false;       
           end
       end
       function finalize(obj) %#ok<INUSL>
           obj.ObjFun.Correspondence.Visible = false;
       end
       function Mstep(obj)
           updateCorrespondence(obj.ObjFun);% updating correspondences
           match(obj.ObjFun.Tmodel,obj.ObjFun.Correspondence,obj.ObjFun.Floating,obj.ObjFun.B);% Weighted Update of Tmodel
           if obj.LevelUpdate, return; end% specialized eval required
           eval(obj.ObjFun);% eval current situation
       end
   end
end % classdef