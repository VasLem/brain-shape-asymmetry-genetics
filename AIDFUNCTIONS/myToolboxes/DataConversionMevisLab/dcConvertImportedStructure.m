function [out] = dcConvertImportedStructure(input,outputpath)
    n = length(input);
    [path,ID] = setupParForProgress(n);
    bk = pwd;
    if ~strcmp(outputpath(end),'/'),outputpath = [outputpath '/'];end
    outputpath = [outputpath 'IMAGES/'];
    out = zeros(8,n);
    parfor i=1:1:n
        tmpout = zeros(8,1);
        tmpout(1) = 1;
        try
           % i = 10;      
           cd(input{i});
           infofile = dir('*.info');
           if isempty(infofile), out(:,i) = tmpout; continue; end 
           A=fileread(infofile.name);
           B=strfind(A,'originalfilename'); %B gives indices of start 'originalfilename' string
           if isempty(B), out(:,i) = tmpout; continue; end
           C=strfind(A,newline); %detects when new line starts
           id=find((C-(B+19))>0); %detect end orignal filename string
           originalname=A(B+19:C(id(1))-1);  %original filename starts then on posistion B + 19 
           %Save to points and cleaned folder
           if ~ischar(originalname), out(:,i) = tmpout; continue; end
           tmpout(2) = 1;
           if ~isempty(dir('initialization.xml'))
              copyfile('initialization.xml',[outputpath '20 XML POSE POINTS/' originalname '.xml']);
              posepoints = readXMLLandmarks('initialization.xml');
              saveLandmarkFile(posepoints,originalname,outputpath);
              tmpout(3) = 1;
           end
           if ~isempty(dir('original.obj'))
              copyfile('original.obj',[outputpath '01 ORIGINAL IMAGES/' originalname '.obj']);
              tmpout(4) = 1;
           end
           if ~isempty(dir('original.bmp'))
              copyfile('original.bmp',[outputpath '01 ORIGINAL IMAGES/' originalname '.bmp']);
              tmpout(5) = 1;
           end
           if ~isempty(dir('original.gif'))
              copyfile('original.gif',[outputpath '01 ORIGINAL IMAGES/' originalname '.gif']);
              tmpout(6) = 1;
           end
           if ~isempty(dir('original.mtl'))
              copyfile('original.mtl',[outputpath '01 ORIGINAL IMAGES/' originalname '.mtl']);
              tmpout(7) = 1;
           end
           if ~isempty(dir('purified.wem'))
              copyfile('purified.wem',[outputpath '30 WEM CLEANED/' originalname '.wem']);
              tmpout(8) = 1;
           end
        catch
            tmpout(1) = 0;
        end
        out(:,i) = tmpout;
        cd(bk);
       parfor_progress;
    end
    closeParForProgress(path,ID);
end

% out is 7 by N files
% 1-> failed or not, 2 -> info, 3 -> xml, 4 ->original.obj, 5->
% orginal.bmp, 6 ->original.gif, 7-> orginal.mtl, 8 -> purified.wem


function saveLandmarkFile(posepoints,originalname,outputpath) %#ok<*INUSL>
    save([outputpath '21 MAT POSE POINTS/' originalname '.mat'],'posepoints');
end