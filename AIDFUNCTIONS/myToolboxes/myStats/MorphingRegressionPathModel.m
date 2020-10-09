function MorphingRegressionPathModel(Model,P,range,step,rec)
         if nargin < 5, rec = false; end
         cl = clone(Model.Average);
         v = viewer(cl);
         cl.Material = 'Dull';
         cl.SingleColor = [0.7 0.7 0.7];
         v.SceneLightVisible = true;
         Y = Model.AvgCoeff;
         MovieFile = [];
         if rec % initialize recording
            [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.avi','Avi files (*.avi)'}, ...
                                              'Save as', 'Movie'); %#ok<NASGU>
            if filename == 0, 
               rec = false; 
            else
                  cd(pathname);
                  MovieFile = avifile(filename);
                  MovieFile.fps = 5;
                  MovieFile.compression = 'none';
               end
         end         
         in = input('Ready? y/n:','s');
         while strcmp(in,'y')
             for i=range(1):step:range(2)
                 DY = P*i;
                 NY = Y+DY';
                 updateShowPC(Model,cl,NY);
                 if rec,MovieFile = recordFrame(MovieFile,v);end
                 pause(0.1);
             end
             pause;
             for i=range(2):-step:range(1)
                 DY = P*i;
                 NY = Y+DY';
                 updateShowPC(Model,cl,NY);
                 if rec,MovieFile = recordFrame(MovieFile,v);end
                 pause(0.1);
             end
                 in = input('Again? y/n:','s');
         end
         if rec, MovieFile = close(MovieFile);end
         delete(v);
end

function movie = recordFrame(movie,f)
                      switch class(f)
                          case 'double'
                             set(f,'Resize','off');
                             im = getframe(f);
                          case 'viewer3DObj'
                             set(f.Figure,'Resize','off');
                             im = getframe(f.Figure);
                      end                      
                      movie = addframe(movie,im);
             end