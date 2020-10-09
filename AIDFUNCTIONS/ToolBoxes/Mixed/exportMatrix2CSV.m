function exportMatrix2CSV(M,filename,path)
         currentpath = pwd;
         if nargin==3, cd(path);end
         fileID = fopen([filename '.txt'],'w');
         [s1,s2]=size(M);
         formatSpec = [];
         for i=1:1:s2-1
             formatSpec = [formatSpec '%2.6e,'];
         end
         formatSpec = [formatSpec '%2.6e\n'];
         for i=1:1:s1
             fprintf(fileID,formatSpec,M(1,:));
         end
         fclose(fileID);
         cd(currentpath);
end