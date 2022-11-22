classdef System < matlab.System

    % == TUNABLE PROPERTIES ==
    properties
        m_x = 40;
        m_y = 3;
        y_max = 1.2;
        g = 9.81;
        c_varphi = 0.12;
        c_x = 0.5;
        k_x = 300;
        k_l = 500;
        l_0 = 0.23;
        k_vx = 3.6;
        k_vl = -1.65;
        k_phi = 0.4;
        
    end

    % == DISCRETE STATES ==
    properties (DiscreteState)
        L
        X
        l
        phi
        varphi
        x
        
    end

    % == METHODS ==
    methods (Access = protected)

        % -- Initialization of the states --
        function resetImpl(obj)
            
            obj.L = 0.0;
            obj.X = 0.0;
            obj.l = 0.0;
            obj.phi = 0.0;
            obj.varphi = 0.0;
            obj.x = 0.0;
            
        end

        % -- Implementation of actual system behavior --
        function [d_L, d_X, d_l, d_phi, d_varphi, d_x,x_out, l_out, phi_out] = stepImpl(obj, L, X, l, phi, varphi, x, v_x, v_l)
            % This is an artifact of the code generation process. The sympy code generation couldn't
            % possibly know that the additional properties are only available as object properties, which is
            % we map them to local variables here, which will then be usable within the automatically
            % generated state equation expressions.
            m_x = obj.m_x;
            m_y = obj.m_y;
            y_max = obj.y_max;
            g = obj.g;
            c_varphi = obj.c_varphi;
            c_x = obj.c_x;
            k_x = obj.k_x;
            k_l = obj.k_l;
            l_0 = obj.l_0;
            k_vx = obj.k_vx;
            k_vl = obj.k_vl;
            k_phi = obj.k_phi;
            

            % Here we update the internal state from exactly the inputs
            obj.L = L;
            obj.X = X;
            obj.l = l;
            obj.phi = phi;
            obj.varphi = varphi;
            obj.x = x;
            

            % First of all we calculate the output using the current states and the given input
            x_out = x;
            l_out = l;
            phi_out = k_phi * rad2deg(varphi) ;
            
            v_x = k_vx * v_x;
            v_l = k_vl * v_l;

            if (x < 0 && X < 0) || (x > 2.5 && X > 0)
                v_x = 0;
            end

            if (l < 0 && L < 0) || (l > 1.3 && L > 0)
                v_l = 0;
            end

            % Then we apply the state equations so that we can output the derived system state, which will
            % then have to be fed into an integrator block
            d_L = L.*k_l.*l.*m_x./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - L.*k_l.*l.*m_y.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + 2*L.*k_l.*l.*m_y./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + L.*k_l.*l_0.*m_x./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - L.*k_l.*l_0.*m_y.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + 2*L.*k_l.*l_0.*m_y./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + X.*k_x.*l.*m_y.*sin(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + X.*k_x.*l_0.*m_y.*sin(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + c_varphi.*m_y.*phi.*sin(varphi).*cos(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + g.*l.*m_y.^2.*sin(varphi).^2.*cos(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + g.*l_0.*m_y.^2.*sin(varphi).^2.*cos(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - k_l.*l.*m_x.*v_l./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + k_l.*l.*m_y.*v_l.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 2*k_l.*l.*m_y.*v_l./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - k_l.*l_0.*m_x.*v_l./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + k_l.*l_0.*m_y.*v_l.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 2*k_l.*l_0.*m_y.*v_l./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - k_x.*l.*m_y.*v_x.*sin(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - k_x.*l_0.*m_y.*v_x.*sin(varphi)./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - l.^2.*m_x.*m_y.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + l.^2.*m_y.^2.*phi.^2.*sin(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + l.^2.*m_y.^2.*phi.^2.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 2*l.^2.*m_y.^2.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 2*l.*l_0.*m_x.*m_y.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + 2*l.*l_0.*m_y.^2.*phi.^2.*sin(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + 2*l.*l_0.*m_y.^2.*phi.^2.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 4*l.*l_0.*m_y.^2.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - l_0.^2.*m_x.*m_y.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + l_0.^2.*m_y.^2.*phi.^2.*sin(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) + l_0.^2.*m_y.^2.*phi.^2.*cos(varphi).^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2) - 2*l_0.^2.*m_y.^2.*phi.^2./(-l.*m_x.*m_y + l.*m_y.^2.*sin(varphi).^2 + l.*m_y.^2.*cos(varphi).^2 - 2*l.*m_y.^2 - l_0.*m_x.*m_y + l_0.*m_y.^2.*sin(varphi).^2 + l_0.*m_y.^2.*cos(varphi).^2 - 2*l_0.*m_y.^2);
            d_X = L.*k_l.*l.*sin(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + L.*k_l.*l_0.*sin(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + X.*k_x.*l./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + X.*k_x.*l_0./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + c_varphi.*phi.*cos(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + g.*l.*m_y.*sin(varphi).*cos(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) + g.*l_0.*m_y.*sin(varphi).*cos(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) - k_l.*l.*v_l.*sin(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) - k_l.*l_0.*v_l.*sin(varphi)./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) - k_x.*l.*v_x./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y) - k_x.*l_0.*v_x./(-l.*m_x + l.*m_y.*sin(varphi).^2 + l.*m_y.*cos(varphi).^2 - 2*l.*m_y - l_0.*m_x + l_0.*m_y.*sin(varphi).^2 + l_0.*m_y.*cos(varphi).^2 - 2*l_0.*m_y);
            d_l = L;
            d_phi = L.*k_l.*l.*m_y.*sin(varphi).*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + L.*k_l.*l_0.*m_y.*sin(varphi).*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 2*L.*l.*m_x.*m_y.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - 2*L.*l.*m_y.^2.*phi.*sin(varphi).^2./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - 2*L.*l.*m_y.^2.*phi.*cos(varphi).^2./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 4*L.*l.*m_y.^2.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 2*L.*l_0.*m_x.*m_y.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - 2*L.*l_0.*m_y.^2.*phi.*sin(varphi).^2./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - 2*L.*l_0.*m_y.^2.*phi.*cos(varphi).^2./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 4*L.*l_0.*m_y.^2.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + X.*k_x.*l.*m_y.*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + X.*k_x.*l_0.*m_y.*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + c_varphi.*m_x.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - c_varphi.*m_y.*phi.*sin(varphi).^2./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 2*c_varphi.*m_y.*phi./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + g.*l.*m_x.*m_y.*sin(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - g.*l.*m_y.^2.*sin(varphi).^3./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 2*g.*l.*m_y.^2.*sin(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + g.*l_0.*m_x.*m_y.*sin(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - g.*l_0.*m_y.^2.*sin(varphi).^3./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) + 2*g.*l_0.*m_y.^2.*sin(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - k_l.*l.*m_y.*v_l.*sin(varphi).*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - k_l.*l_0.*m_y.*v_l.*sin(varphi).*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - k_x.*l.*m_y.*v_x.*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2) - k_x.*l_0.*m_y.*v_x.*cos(varphi)./(-l.^2.*m_x.*m_y + l.^2.*m_y.^2.*sin(varphi).^2 + l.^2.*m_y.^2.*cos(varphi).^2 - 2*l.^2.*m_y.^2 - 2*l.*l_0.*m_x.*m_y + 2*l.*l_0.*m_y.^2.*sin(varphi).^2 + 2*l.*l_0.*m_y.^2.*cos(varphi).^2 - 4*l.*l_0.*m_y.^2 - l_0.^2.*m_x.*m_y + l_0.^2.*m_y.^2.*sin(varphi).^2 + l_0.^2.*m_y.^2.*cos(varphi).^2 - 2*l_0.^2.*m_y.^2);
            d_varphi = phi;
            d_x = X;
            

        end

    end

end