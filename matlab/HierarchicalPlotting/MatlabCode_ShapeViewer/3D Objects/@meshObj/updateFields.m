function updateFields(obj,action,varargin)
% updating fields according to action, typically called when changing faces
% and vertices
         switch action
             case 'Reduce' % varargin{1} = index for values to keep
                 keepindex = varargin{1};
                 if isempty(keepindex), updateFields(obj,'Clear'); return; end
                 %if ~isempty(obj.Gradient), obj.Gradient = obj.Gradient(:,keepindex);end
                 try 
                    if ~isempty(obj.IndexedColor), obj.IndexedColor = obj.IndexedColor(:,keepindex);end
                 catch
                 end
                 if isempty(obj.TextureMap)||isempty(obj.UV),
                    if ~isempty(obj.TextureColor), obj.TextureColor = obj.TextureColor(:,keepindex);end
                 end
%                  if ~isempty(obj.Distance)
%                     if ~length(keepindex)==length(obj.Distance),obj.Distance = obj.Distance(:,keepindex);end
%                  end
                 try
                     obj.Distance = obj.Distance(:,keepindex);
                 catch %#ok<CTCH>
                 end
                 if ~isempty(obj.UV), obj.UV = obj.UV(:,keepindex);end
                 updateChildren(obj,'Reduce',varargin{:});                
             case 'Clear'
                 delete(obj.ph);
                 %obj.Gradient = [];
                 obj.TextureColor = [];
                 obj.IndexedColor = [];
                 obj.Distance = [];
                 obj.UV = [];
                 updateChildren(obj,'Delete');
             case 'Append'
                 %if ~isempty(obj.Gradient), obj.Gradient = [obj.Gradient varargin{1}.Gradient]; end
                 if ~isempty(obj.IndexedColor), obj.IndexedColor = [obj.IndexedColor varargin{1}.IndexedColor]; end
                 if ~isempty(obj.TextureColor), obj.TextureColor = [obj.TextureColor varargin{1}.TextureColor]; end
                 if ~isempty(obj.UV), obj.UV = [obj.UV varargin{1}.UV]; end
                 if ~isempty(obj.Distance), obj.Distance = [obj.Distance varargin{1}.Distance]; end
                 updateChildren(obj,'Append');
             otherwise
                     return;
         end
end