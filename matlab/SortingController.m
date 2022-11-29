classdef SortingController < matlab.System
    % untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        positions = [0.3 0.5 0.7 0.9 1.1];
        masses = [1 7 3 4 6];
        ordering = [5 2 1 3 4];
        pickup_length = 1.0;
        dropoff_length = 0.1;
        dropoff_position = 2.0

        box_height = 0.2;
        
        % The speed and the time for the raising of the crane at pick up 
        % or drop off

        crane_ascend_speed = 0.5;
        crane_ascend_time = 1.0;

        % The speed and the time for the raisung if the crane at pick up 
        % or drop off

        crane_descend_speed = 0.05;
        crane_descend_time = 0.5;

        % This is the time where all the measurements have to be within 
        % the tolerance margins to declare the final state to be reached
        
        t_tol = 1.0;

        % these are the tolerances which are allowed for the conditional 
        % checks of whether the destination has been reached

        x_tol = 0.05;
        l_tol = 0.05;
        phi_tol = 1.0;

        % This is a boolean flag which needs to be set when doing
        % simulations not working with the real model. The reason for this
        % is that the output measurments of the simulated model have some
        % weird jitter which break the condition to check if the model 
        % has arrived in the destination position. If this flag is active 
        % we use a reduced condition which works in simulation as well
        simulation_fix_condition = 1;
    end

    properties (DiscreteState)
        % This is the integer representation of the current program state 
        % (basically which phase we are currently in)
        % 0 - "Go to location"
        % 1 - "lower crane"
        % 2 - "raise crane"
        state

        % this integer state variable identifies if we are currently 
        % picking up or dropping off
        % 0 - "pick up"
        % 1 - "drop off"
        transport_state
           
        % This is a special flag which has to be set whenenver a state 
        % transition is being done to signal that in that very same 
        % iteration no other states should be activated
        skip

        % This variable saves the integer index of the how many-th box 
        % we are currently processing. For example if we are just at the 
        % beginning of the program then this will be 1, because we process 
        % the first box, later on we increment it as we go through all the 
        % 5 boxes.
        index

        % in the end we want to stack the boxes in a certain order which 
        % we are able to determine ourselves. So this index will save the 
        % actual index which doesn't just monotonically increases but 
        % instead represents this order
        ordered_index

        % This is the boolean state of whether or not the current 
        % position measurements are within the required margin at the 
        % destination position which is needed for pick up and drop off.
        within_margin

        dest_length

        % This is a buffer value where we can save a timestamp. For various
        % program states we need a duration since the beginning of that 
        % phase and saving the time at the beginning of the phase into this
        % variable is how we achieve that
        start_time

        % This variable simply stores the external clock time internally
        time

        x_internal
        l_internal
        v_x_internal
        v_l_internal
        mag_internal
        control_internal
        m_internal
    end

    % Pre-computed constants
    properties (Access = private)

    end

    methods (Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end
        
        %% == INITIALIZATION ==
        % In this method we have to assign all the initial values to the 
        % internal states
        function resetImpl(obj)
            obj.state = 0.0;
            obj.transport_state = 0.0;
            
            obj.index = 1;
            obj.ordered_index = obj.ordering(obj.index);
            obj.dest_length = obj.pickup_length;
            obj.start_time = 0.0;
            obj.within_margin = 0.0;
            
            obj.time = 0.0;

            obj.x_internal = 0.0;
            obj.l_internal = 0.0;
            obj.v_x_internal = 0.0;
            obj.v_l_internal = 0.0;
            obj.mag_internal = 0.0;
            obj.control_internal = 0.0;
            obj.m_internal = 0.0;
        end

        %% == STEP ==
        % This method is called in every iteration of the simulink 
        % simulation. If the simulation step size is for example T=0.01 
        % then this method is called every 0.01 seconds.
        function [x, l, v_x, v_l, mag, control, m, state] = stepImpl(obj, time, x_mess, l_mess, phi_mess, reset)
            
            obj.time = time;
            duration = time - obj.start_time;
            obj.skip = 0;
            
            if obj.index < 6
                obj.ordered_index = obj.ordering(obj.index);
            end

            % termination condition: If the last box was processed simply 
            % stop.
            if obj.index > size(obj.positions)
                obj.state = 10.0;
            end

            % If we give the external reset signal, then we reset the 
            % internal sorting state most importantly, but also as long 
            % as the signal is active the crane and the gantry slowly move
            % back to the original position.
            if reset > 0.9 
                % reset the internal states
                obj.state = 0;
                obj.transport_state = 0;
                obj.index = 1;
                obj.skip = 1;

                obj.control_internal = 0;
                obj.v_x_internal = 0;
                obj.v_l_internal = 0;

                if x_mess > 0.05
                    obj.v_x_internal = -0.05;
                end

                if l_mess > 0.05
                    obj.v_l_internal = -0.05;
                end
                
            end

            % State 0 -- Go to location
            % In this state the control is engaged and we use the position 
            % control to go to the pick up location.
            if obj.state == 0.0 && ~(obj.skip == 1)
                obj.control_internal = 1.0;

                if (obj.transport_state == 0.0)
                    x_dest = obj.positions(obj.ordered_index);
                    l_dest = obj.pickup_length;
                else
                    x_dest = obj.dropoff_position;
                    l_dest = obj.dest_length;
                end

                obj.x_internal = x_dest;
                obj.l_internal = l_dest;
                
                % Now we need to somehow make sure that we are actually 
                % at the destination and that the position is stable 
                % as well. We do that by checking if we can stay within 
                % A boundry for a certain amount of time.
                x_dist = abs(x_dest - x_mess);
                l_dist = abs(l_dest - l_mess);
                phi_dist = abs(phi_mess);

                % The lower condition does not work in simulation for some 
                % weird reason...
                if obj.simulation_fix_condition
                    condition = (phi_dist < obj.phi_tol);
                else
                    condition = (x_dist < obj.x_tol) && (l_dist < obj.l_tol) && (phi_dist < obj.phi_tol);
                end

                if condition
                    if (obj.within_margin == 0.0)
                        obj.within_margin = 1.0;
                        obj.start_time = time;
                    else
                        if (duration > obj.t_tol)
                            % STATE TRANSITION -- Lower crane                            
                            obj.state = 1.0;
                            obj.start_time = time;
                            obj.skip = 1;
                        end
                    end
                else
                    obj.within_margin = 0;
                end
            end
            

            % State 1 -- Lower crane
            % In this state we lower the crane with a fixed speed for a 
            % certain amount of time. This could be either to pick up 
            % an item or to drop one off.
            if obj.state == 1.0 && ~(obj.skip == 1)
                obj.control_internal = 0.0;
                obj.v_x_internal = 0.0;
                obj.v_l_internal = obj.crane_descend_speed;

                if duration > obj.crane_descend_time 
                   % TRANSITION 2 - Ascend crane
                   obj.state = 2.0;
                   obj.start_time = obj.time;


                   % only if we are currently dropping off we disengage 
                   % the magnet at the bottom otherwise we engage
                   if obj.transport_state == 1.0
                        obj.mag_internal = 0.0;
                        obj.m_internal = 0.0;
                        obj.transport_state = 0.0;
                        obj.dest_length = obj.dest_length - obj.box_height;
                   else
                        obj.mag_internal = 1.0;
                        obj.m_internal = obj.masses(obj.ordered_index);
                        obj.transport_state = 1.0;
                   end
                end
            end


            % State 2 -- Raise crane
            % In this state we raise the crane with a fixed speed and a 
            % fixed amount of time. This could either at pick up or drop 
            % off to avoid collision with other boxes.
            if obj.state == 2.0 && ~ (obj.skip == 1)
                obj.control_internal = 0.0;
                obj.v_x_internal = 0.0;
                obj.v_l_internal = -obj.crane_ascend_speed;
            
                if duration > obj.crane_ascend_time
                    % TRANSITION 0 - move to position
                    obj.state = 0.0;
                    obj.start_time = time;

                    % only if we are currently in "pick up" mode, which 
                    % means this is the ascend directly following the 
                    % a drop off at the destination, we increase the 
                    % index by one which means that we now process the next
                    % box in the array.
                    if obj.transport_state == 0.0
                        obj.index = obj.index + 1;
                        
                    end
                end
            end

            % STATE 10 - terminal state where we do not do anything
            if obj.state == 10.0
                obj.control_internal = 0.0;
                obj.v_l_internal = 0.0;
                obj.v_x_internal = 0.0;
            end

            % Up to this point we have only ever modified the internal 
            % state equivalents of the output values. Now we need to
            % actually set the output variables to the same values as the 
            % these internal representations.
            x = obj.x_internal;
            l = obj.l_internal;
            v_x = obj.v_x_internal;
            v_l = obj.v_l_internal;
            mag = obj.mag_internal;
            control = obj.control_internal;
            m = obj.m_internal;
            state = obj.state;
        end
        
        %% == UTILITY METHODS ==

        function stateTransition(obj, state)
            obj.state = state;
            obj.start_time = obj.time;
        end
        
    end
end
