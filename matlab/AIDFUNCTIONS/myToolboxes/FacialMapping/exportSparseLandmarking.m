function exportSparseLandmarking(DATA,LandmarkNames,filename)
         headers = {'ID' 'LM' 'Name' 'Side' 'X1' 'Y1' 'Z1' 'X2' 'Y2' 'Z2' 'X3' 'Y3' 'Z3'};
         nheaders = length(headers);
         N = length(DATA);
         nLM = DATA{1}.LM{1,2}.nVertices;
         %  AIM IS TO BUILD 4 DATA CELL ARRAYS TO EXPORT
         for o=1:1:2
             %o=1;
             celldata = headers;
             for i=1:1:N
                 %i=1;
                 tmp = DATA{i};
                 tmpcell = cell(nLM,nheaders);
                 tmpcell(:,1) = {num2str(tmp.ID)};
                 tmpcell(:,2) = num2cell(1:nLM)';
                 tmpcell(:,3:4) = LandmarkNames;
                 tmpcell(:,5:7) = num2cell(tmp.LM{o,1}.Vertices);
                 tmpcell(:,8:10) = num2cell(tmp.LM{o,2}.Vertices);
                 tmpcell(:,11:13) = num2cell(tmp.LM{o,3}.Vertices);
                 celldata = [celldata; tmpcell];
             end
             xlswrite(filename,celldata,['Observer ' num2str(o)],'A1');
         end
         % Now exporting the average of the indicators
         celldata = headers;
         for i=1:1:N
             %i=1;
             tmp = DATA{i};
             tmpcell = cell(nLM,nheaders);
             tmpcell(:,1) = {num2str(tmp.ID)};
             tmpcell(:,2) = num2cell(1:nLM)';
             tmpcell(:,3:4) = LandmarkNames;
             tmpcell(:,5:7) = num2cell(tmp.AvgLM{1}.Vertices);
             tmpcell(:,8:10) = num2cell(tmp.AvgLM{2}.Vertices);
             tmpcell(:,11:13) = num2cell(tmp.AvgLM{3}.Vertices);
             celldata = [celldata; tmpcell];
         end
         xlswrite(filename,celldata,'Observer Averages','A1');
         % Now automatic indications
         % Now exporting the average of the indicators
         celldata = headers;
         for i=1:1:N
             %i=1;
             tmp = DATA{i};
             tmpcell = cell(nLM,nheaders);
             tmpcell(:,1) = {num2str(tmp.ID)};
             tmpcell(:,2) = num2cell(1:nLM)';
             tmpcell(:,3:4) = LandmarkNames;
             tmpcell(:,5:7) = num2cell(tmp.AvgLMOnFace{1}.Vertices);
             tmpcell(:,8:10) = num2cell(tmp.AvgLMOnFace{2}.Vertices);
             tmpcell(:,11:13) = num2cell(tmp.AvgLMOnFace{3}.Vertices);
             celldata = [celldata; tmpcell];
         end
         xlswrite(filename,celldata,'Automatic indications','A1');
end