function varargout = importGetFile
         [filename, pathname,filterindex] = uigetfile({'*.obj','Obj Wavefront';...
                                                       '*.obj','Obj 3DMD(2Pod)';...
                                                        },'MultiSelect','on');
         if isequal([filename,pathname],[0,0]),varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4} = [];return; end
         if ~iscell(filename)
            varargout{1} = 1;
            varargout{2} = {filename};
         else
            varargout{1} = size(filename,2);
            varargout{2} = filename;
         end
         varargout{3} = pathname;
         switch filterindex
             case 1
                 varargout{4} = 'Wavefront';
             case 2
                 varargout{4} = '3DMD 2pod';
             otherwise
         end
%          cd(pathname);
end
