function out = myCSVRead(filename)    
      fid = fopen(filename);
      counter = 0;
      while 1
          %disp(counter);
          tline = fgetl(fid);
          if ~ischar(tline), break, end
          counter = counter+1;
          commas = strfind(tline,',');
          commas = [0 commas];
          if counter==1
             nrCol = length(commas); 
             out = cell(1000,nrCol); 
          end
          if mod(counter,1000)==0
             out = [out; cell(1000,nrCol)];
          end
          forcell = cell(1,nrCol);
          for j=1:1:length(commas)
               if j<length(commas)
                 forcell{j} = tline(commas(j)+1:commas(j+1)-1);
               else
                 forcell{j} = tline(commas(j)+1:end);
               end
          end
          out(counter,:) = forcell;
      end
      fclose(fid);
      out = out(1:counter,:);
end