function generateSpokesPlot(fig,X,mu,s)
         nrX = length(X);
         %deg = floor(360/nrX);
         Z = (X-mu)./s;
         xplace = 1:nrX;
         bar(xplace,Z',1)
             


end