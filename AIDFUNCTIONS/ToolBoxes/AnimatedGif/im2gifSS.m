 fn = ls('model*vy-45.bmp');
 %filename template for bmp files
 [a,b] = size(fn);

 
 
for f = 1:a %Only one so far
    fname = fn(f,:)
    ind = findstr(fname,'vy')+1;
    rootname = fname(1:ind);
    gifname = [rootname,'.gif']
    imcnt =  1;    
    
    i = [1,2,3,4,5,6,7,7,6,5,4,3,2,1];
    r = rand(7);  % or however many bmp views there are
    i2 = [i(r:length(i)),i(1:r-1)] %start in a random place
    
    for i = i2  % or [1,2,3,4,5,6,7,7,6,5,4,3,2,1]
        imname = [rootname, num2str(60-i*15), '.bmp'];
       
        inimg = imread(imname);
         [img,imap] = rgb2ind(inimg,256);
         %delayTime should allow you to adjust the speed.  Currently 8 fps
         if imcnt == 1
             imwrite(img, imap, gifname,'delayTime',0.250,'loopCount',Inf,'WriteMode','overwrite');
         else
             imwrite(img, imap, gifname,'delayTime',0.250,'WriteMode','append');
         end
         imcnt = imcnt+1;
    end
    
end
