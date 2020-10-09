function RF = projectStudies(pc,lc,F,S)
          RProj=pc*diag(lc.^(-1/2));
          %[a,b] = size(F);
          [a,b] = size(RProj);
          minF=min(F);
          minF = minF(1:b);
          maxF=max(F);
          maxF = maxF(1:b);
          RF = zeros(a,b,size(S,3));
          for k=1:1:size(S,3)
              RF(:,:,k) = S(:,:,k)*RProj;
              minF=min([minF;RF(:,:,k)]);
              maxF=max([maxF;RF(:,:,k)]);
          end
          figure;hold on;col = colormap('lines');
%           ind = randperm(64);
%           r = ind(1:size(S,3));
          r = randsample(64,size(S,3),true);
          colors = col(r,:);
          names = cell(1,a);
          for i=1:1:a
              names{i} = num2str(i);
          end
          
%     figure;clf;
%     ax1=1;ax2=2;
%     plotxyrange(F,ax1,ax2,...
%      'Compromise as barycenter of Studies',char(names),...
%     [minF(ax1),maxF(ax1),minF(ax2),maxF(ax2)]);

   plot3(F(:,1),F(:,2),F(:,3),'m.','MarkerSize',30);hold on;
%    for k=1:size(S,3);
%        %k=2;
%        G = squeeze(RF(:,:,k));
%        plot3(G(:,1),G(:,2),G(:,3),'Color',colors(k,:),'LineStyle','none','Marker','.');
%        for i=1:a;
%            % i=1;
%            x=[F(i,ax1),G(i,ax1)];
%            y=[F(i,ax2),G(i,ax2)];
%            plot(x,y,'-',...
%                  'Color',colors(k,:));
%        end
%    end         
end