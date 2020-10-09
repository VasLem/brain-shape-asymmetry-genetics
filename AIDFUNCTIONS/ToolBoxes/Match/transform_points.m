function [new_points] = transform_points(points, T)
        
        if size(points,2)==3
            points = points';
        end

        %points
        tmp = [points; ones(1,size(points,2))];
        tmp = T * tmp;
        new_points = tmp(1:3,:);
        