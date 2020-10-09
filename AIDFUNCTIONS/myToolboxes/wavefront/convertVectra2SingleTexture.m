function [map,convert] = convertVectra2SingleTexture(maps)
         w = zeros(1,3);
         h = zeros(1,3);
         for i=1:3
             [w(i),h(i),~] = size(maps{i+1});
         end
         %make all images same w and h
         W = max(w);
         H = max(h);
         IM =0*ones(W,H,3,'uint8');
         
         convert = cell(1,4);
         convert{1} = [0 0;0 0];
         im = cell(1,3);
         for i=1:1:3
            convert{i+1} = [0 1;0 1];
            tmpw = w(i);
            tmph = h(i);
            
            convert{i+1}(1,2) = tmph/H;
            %convert{i+1}(2,2) = tmpw/W;
            convert{i+1}(2,1) = 1-tmpw/W;

            im{i} = IM;
            im{i}(1:tmpw,1:tmph,:) = maps{i+1};
            %figure;imshow(im{i})
         end
         
         convert{1+1} = convert{1+1}/2;
         
         convert{2+1}(1,:) = convert{2+1}(1,:)/2;
         convert{2+1}(2,:) = convert{2+1}(2,:)/2+0.5;
         
         convert{3+1}(1,:) = convert{3+1}(1,:)/2+0.5;
         convert{3+1}(2,:) = convert{3+1}(2,:)/2;
         
         %W = matrix ROWS = top row in convert;
         %H = matrix COLLOMNS = bottom row in convert;
         
%          map = 0*ones(2*W,2*H,3,'uint8');
%          map(1:W,1:H,:) = im{1};
%          map(W+1:2*W,1:H,:) = im{2};
%          map(1:W,H+1:2*H,:) = im{3};
         map = 0*ones(2*W,2*H,3,'uint8');
         map(W+1:2*W,1:H,:) = im{1};
         map(1:W,1:H,:) = im{2};
         map(W+1:2*W,H+1:2*H,:) = im{3};
end