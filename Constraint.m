classdef Constraint
    properties
        center
        size
        duration
    end
    methods
        function this = Constraint(center, size, duration)
            this.center = center; % in m?
            this.size = size;
            this.duration = duration; % in ms?
        end
    end
end

