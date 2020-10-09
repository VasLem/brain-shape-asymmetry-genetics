function [TransformedVertices, T] = my_TPS_nonrigid_3D(StartVertices, EndVertices,FloatingVertices,transformSmooth)
         

         type = 'thinplate';
         transformConstant = 3;
         %transformSmooth = 0;

         % FIT
         
         %X-DISP                       
         T.X = rbfcreate(StartVertices',EndVertices(:,1)','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);
         %Y-DISP                       
         T.Y = rbfcreate(StartVertices',EndVertices(:,2)','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);
         %Z-DISP                       
         T.Z = rbfcreate(StartVertices',EndVertices(:,3)','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);            
         
         % EVAL
         X = rbfinterp(FloatingVertices',T.X);
         Y = rbfinterp(FloatingVertices',T.Y);
         Z = rbfinterp(FloatingVertices',T.Z);                                                      
         TransformedVertices = [X(:) Y(:) Z(:)];                       
                                                      
end