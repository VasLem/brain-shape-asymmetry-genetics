function str = myNum2Str(num,padding)
         str = num2str(num);
         n = length(str);
         if n>=padding, return; end% no zeros to add
         add = [];
         for i=1:padding-n
             add = [add '0']; %#ok<AGROW>
         end
         str = [add str];
         return;
end