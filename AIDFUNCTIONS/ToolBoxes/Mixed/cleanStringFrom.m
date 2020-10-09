function str = cleanStringFrom(str,from)
      % underscore
          index = (1:length(str));
          bad = strfind(str,from);
          str = str(setdiff(index,bad));
end