function out = replaceColor(obj,varargin)
         if nargout == 1;
             obj = clone(obj);
             out = obj;
         end
         Input = find(strcmp(varargin,'Old'));
         if ~isempty(Input)
             Old = varargin{Input+1};
         else
             Old = impixel(obj);
         end
         Input = find(strcmp(varargin,'New'));
         if ~isempty(Input)
             New = varargin{Input+1};
         else
             New = colorui;
             if ~(length(color) == 1), return; end;
         end
         Input = find(strcmp(varargin,'Distance'));
         if ~isempty(Input)
             Dist = varargin{Input+1};
         else
             Dist = 0;
         end
         % Red Channel
         index = find((obj.Red>=Old(1)-Dist)&(obj.Red<=Old(1)+Dist)&...
                      (obj.Green>=Old(2)-Dist)&(obj.Green<=Old(2)+Dist)&...
                      (obj.Blue>=Old(3)-Dist)&(obj.Blue<=Old(3)+Dist));
         obj.Red(index) = New(1);
         obj.Green(index) = New(2);
         obj.Blue(index) = New(3);
end