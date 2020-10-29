function [usemtl,mtl,maps] = parseMTL(filename,usePrefix)
        if nargin<2, usePrefix = []; end
        usemtl = false;maps = [];mtl = [];
        MTLName = [filename(1:end-4) '.mtl'];
        index = strfind(MTLName,'/');
        if isempty(index)
           path = pwd;
        else
           path = MTLName(1:index(end));
        end
        tmp = dir(MTLName);
        if isempty(tmp), return;end
        usemtl = true;
        fid = fopen(MTLName,'r');
        MTL = textscan(fid,'%[^\n\r]');MTL = MTL{1};
        fclose(fid);
        nlines = length(MTL);
        % searching for mtl definition(s)
        mtl = {};
        newmtlindex = [];
        for i=1:1:nlines
            str = MTL{i};
            if ~contains(str,'newmtl'), continue; end
            index = strfind(str,' ');
            mtl{end+1} = str(index+1:end); %#ok<*AGROW>
            newmtlindex = [newmtlindex i];
        end
        nnewmtl = length(mtl);
        % searching for map_Kd definitions
        map_Kd = cell(1,nnewmtl);
        for i=1:1:nnewmtl
            startindex = newmtlindex(i);
            if i<nnewmtl
               endindex = newmtlindex(i+1);
            else
               endindex = nlines;
            end
            for j=startindex:endindex
                str = MTL{j};
                if ~contains(str,'map_Kd'), continue; end
                index = strfind(str,' ');
                map_Kd{i} = str(index+1:end);
            end
        end
        nMaps = length(map_Kd);
        maps = cell(1,nMaps);
        for i=1:1:nMaps    
            if ~isempty(map_Kd{i})
                if isempty(usePrefix), maps{i} = imread([path map_Kd{i}]); continue; end
                if strcmp(usePrefix, map_Kd{i}(1:length(usePrefix)))
                    maps{i} = imread([path map_Kd{i}]);
                else
                    maps{i} = imread([path usePrefix map_Kd{i}]);
                end
            end
        end     
end