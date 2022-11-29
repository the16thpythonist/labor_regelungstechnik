classdef StateTest < matlab.System

    % == TUNABLE PROPERTIES ==
    properties
        m = 1;
        C = -1;
        
    end

    % == DISCRETE STATES ==
    properties (DiscreteState)
        x;
        z;
        
    end

    % == METHODS ==
    methods (Access = protected)

        % -- Initialization of the states --
        function resetImpl(obj)
            
            obj.x = 0;
            obj.z = 0;
            
        end

        % -- Implementation of actual system behavior --
        function [d_x, d_z,y1, y2] = stepImpl(obj, x, z, u1, u2)
            % This is an artifact of the code generation process. The sympy code generation couldn't
            % possibly know that the additional properties are only available as object properties, which is
            % we map them to local variables here, which will then be usable within the automatically
            % generated state equation expressions.
            m = obj.m;
            C = obj.C;
            
            % First of all we calculate the output using the current states and the given input
            y1 = x ;
            y2 = z ;
            
            % Here we update the internal state from exactly the inputs
            obj.x = x;
            obj.z = z;

            % Then we apply the state equations so that we can output the derived system state, which will
            % then have to be fed into an integrator block
            d_x = m * z;
            d_z = C * x + C * z + u1;
            %d_x = 0;
            %d_z = 0;


        end

    end

end
