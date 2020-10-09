function deletePatch(obj)
         if ~(isempty(obj.ph)) && ishandle(obj.ph), delete(obj.ph); end
end