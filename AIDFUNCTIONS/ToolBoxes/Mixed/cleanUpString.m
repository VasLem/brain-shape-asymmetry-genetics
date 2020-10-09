function out = cleanUpString(in)
      out = in;
      % underscore
          index = (1:length(out));
          bad = strfind(out,'_');
          out = out(setdiff(index,bad));
      % space
          index = (1:length(out));
          bad = strfind(out,' ');
          out = out(setdiff(index,bad));
      % underscore
          index = (1:length(out));
          bad = strfind(out,'.');
          out = out(setdiff(index,bad)); 
      % \
          index = (1:length(out));
          bad = strfind(out,'\');
          out = out(setdiff(index,bad));
      % /
          index = (1:length(out));
          bad = strfind(out,'/');
          out = out(setdiff(index,bad));    
end