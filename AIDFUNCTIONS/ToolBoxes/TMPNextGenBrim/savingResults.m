function savingResults(v)
          %save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end