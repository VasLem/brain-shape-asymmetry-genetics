function captureFrame(obj)
%       f = getframe(obj.Figure);
%       if isempty(obj.MovieFile), return; end
%       obj.MovieFile = addframe(obj.MovieFile,f);
        %figure(obj.Figure);
        %drawnow;
        obj.Status = 'Recording';
        im = captureImage(obj);
%         f = getframe(obj.Figure);
%         im = frame2im(f);
        f = im2frame(im,[]);
        if isempty(obj.MovieFile), return; end
        obj.MovieFile = addframe(obj.MovieFile,f);
%         obj.FrameCounter = obj.FrameCounter +1;        
%         imwrite(im,[num2str(obj.FrameCounter) '.bmp'],'bmp');
        obj.Status = 'Ready';
end