function varargout = normaliseScalePose(obj,refPoseLM,varargin)
         if ~strcmp(class(refPoseLM),'LMObj'), error('No poseLM object'); end        
         if isempty(obj.PoseLM)
            button = questdlg('Would you like to indicate them?','No PoselM!');
            switch lower(button)
                case 'yes'
                    indicatePoseLM(obj);
                otherwise
                    varargout{1} = [];
                    return;
            end
         end
         %Input = strcmp(
         if nargout == 1
             obj = clone(obj);
             obj.Visible = false;
             varargout{1} = obj;
         end
         T = scaledRigidTM;
         match(T,refPoseLM.Vertices,obj.PoseLM.Vertices);
         transform(T,obj);
         delete(T);
end