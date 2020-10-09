function out = chainRule(obj,in,p) %#ok<INUSD>
         if isField(obj,obj.ActiveField)
            out = obj.Scale*obj.Rotation*TM.getPoints(in);
         else
            out = in;% dummy
         end
end