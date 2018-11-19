function printFootSteps(constraints, actTime, steps)
%% function used to print the footsteps
 % in ms?
    
    for i=1:size(constraints,2)
%         constraintsX(i).center = center; % in m?
%         constraintsX(i).size = size;
        constraints(i).duration = duration;
%         if (i*20) < actTime
        x1 = (constraints(i).center(1) - constraints(i).size(1));
        x2 = (constraints(i).center(1) + constraints(i).size(1));
        y1 = (constraints(i).center(2) - constraints(i).size(2));
        y2 = (constraints(i).center(2) + constraints(i).size(2));
        line([x1,x2],[y1,y1]);
        line([x1,x2],[y2,y2]);
        line([x1,x1],[y1,y2]);
        line([x2,x2],[y1,y2]);
%         end
    end
end

