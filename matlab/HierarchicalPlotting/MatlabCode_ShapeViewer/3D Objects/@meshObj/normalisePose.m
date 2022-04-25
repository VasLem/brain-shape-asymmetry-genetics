function varargout = normalisePose(obj,refPoseLM,varargin)
         switch class(refPoseLM)
             case 'LMObj'
                 % ok continue
             case 'meshObj'
                 if isempty(refPoseLM.PoseLM), indicatePoseLM(refPoseLM), end
                 refPoseLM = refPoseLM.PoseLM;
             otherwise
                 error('wrong reference lm');
         end
         %if ~strcmp(class(refPoseLM),'LMObj'), error('No poseLM object'); end        
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
         T = rigidTM;
         match(T,refPoseLM.Vertices,obj.PoseLM.Vertices);
         transform(T,obj);
         delete(T);
end