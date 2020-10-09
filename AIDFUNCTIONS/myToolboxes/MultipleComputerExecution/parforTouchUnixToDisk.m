function parforTouchUnixToDisk(filename)
         unix(['touch ' filename]);
         %save(filename,'in');
end