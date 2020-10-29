function [floatingFeatures, T] = my_TPS_nonrigid_transformation(floatingFeatures, correspondingFeatures,floatingFaces,floatingFlags, inlierWeights,transformSigma)
         

         index = find(inlierWeights>0.5);
         type = 'linear';
         transformConstant = 3;
         transformSmooth = 0;
         
         
         
         % FIT
         
         %X-DISP                       
         T.X = rbfcreate(double(floatingFeatures(index,1:3))',double(correspondingFeatures(index,1))','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);
         %Y-DISP                       
         T.Y = rbfcreate(double(floatingFeatures(index,1:3))',double(correspondingFeatures(index,2))','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);
         %Z-DISP                       
         T.Z = rbfcreate(double(floatingFeatures(index,1:3))',double(correspondingFeatures(index,3))','RBFFunction', type,... 
                         'Stats', 'on','RBFConstant',transformConstant,'RBFSmooth',transformSmooth);            
         
         % EVAL
         X = rbfinterp(floatingFeatures(:,1:3)',T.X);
         Y = rbfinterp(floatingFeatures(:,1:3)',T.Y);
         Z = rbfinterp(floatingFeatures(:,1:3)',T.Z);                                                      
         floatingFeatures(:,1:3) = single([X(:) Y(:) Z(:)]);                       
                                                            
end