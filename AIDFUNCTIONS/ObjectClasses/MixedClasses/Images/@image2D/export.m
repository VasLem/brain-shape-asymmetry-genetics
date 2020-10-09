function export(obj,Field)
          if nargin == 1
             im = obj.Image;
          else
             im = obj.(Field);
          end
          [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.png','Images (*.png)'; ...
                                              '*.bmp','Images (*.bmp)'; ...
                                              '*.jpg','Images (*.jpg)'; ...
                                              '*.gif','Images (*.gif)'},...
                                              'Save as', 'Image');
          if filename == 0; return; end
          cd(pathname);
          switch filterindex
              case 1
                  imwrite(im,filename,'png');
              case 2
                  imwrite(im,filename,'bmp');
              case 3
                  imwrite(im,filename,'jpeg');
              case 4
                  imwrite(im,filename,'gif');
          end
end