function dcPrepareOutputPath(path)
         warning off;
         bk = pwd;
         cd(path);
         mkdir('IMAGES');
         mkdir('RAW');
         cd('IMAGES');
         mkdir('01 ORIGINAL IMAGES');
         mkdir('11 MAT IMAGES');
         mkdir('20 XML POSE POINTS');
         mkdir('21 MAT POSE POINTS');
         mkdir('30 WEM CLEANED');
         mkdir('31 OBJ CLEANED');
         mkdir('32 MAT CLEANED');
         mkdir('41 MAT MAPPED');
         mkdir('42 MAT REFLECTED MAPPED');
         mkdir('51 OBJ MAPPED');
         mkdir('52 OBJ REFLECTED MAPPED');
         cd(bk);
         warning on;
end