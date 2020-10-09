function ViewerAction(obj,what,varargin)
         switch what
             case 'LM Updated'
                 updateTDTable(obj);
             otherwise
         end
end