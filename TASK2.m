classdef TASK2
    properties
        interval; %initial interval
        x_val; %values of the x
        h; %initial value
    end
    methods
        function obj = TASK2()
            obj.interval = [0, 20];
            obj.x_val = [0.4 , 0.3];
            obj.h = 0.8;
        end
        function TASK2aRK4(obj)
            figure(1);
            title("Runge-Kutta 4th order, constant step");
            
            subplot(2,2,1);
            hold on;
            grid on;
            title("x1(t)");
            xlim([0, 20]);
            
            subplot(2,2,2);
            hold on;
            grid on;
            title("x2(t)");
            xlim([0, 20]);
            
            subplot(2,2,3);
            hold on;
            grid on;
            title("x1(t) zoom");
            xlim([12.9, 13.6]);
            ylim([0.7, 0.76]);
            
            subplot(2,2,4);
            hold on;
            grid on;
            title("x2(t) zoom");
            xlim([11.3, 11.9]);
            ylim([0.7, 0.76]);

            legend('show');
            for i = 1:1:10
                [x, t] = Runge_Kutta_Constant(obj, obj.h);
                if obj.h == 0.4
                    figure(5); 
                    title("x1(t) and x2(t) larger step size 0.4");
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'step(0.4): x1(t)');
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'step(0.4): x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    title("x2(x1) step 0.4 and 0.2")
                    plot(x(1, :), x(2, :), 'DisplayName', 'step(0.4): x2(x1)');
                    legend('show');
                    hold on;
                end
                if obj.h == 0.2
                    [t_ode, x_ode] = ode45(@method, obj.interval, obj.x_val);
                    figure(4);
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'x1(t)');
                    title("x1(t) and x2(t) optimal step 0.2 graph + ode");
                    hold on;
                    plot(t_ode, x_ode(:, 1), 'DisplayName', 'ode: x1(t)');
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'x2(t)');
                    hold on;
                    plot(t_ode, x_ode(:, 1), 'DisplayName', 'ode: x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    plot(x(1, :), x(2, :), 'DisplayName', 'x2(x1)');
                    hold on;
                    plot(x_ode(:, 1), x_ode(:, 2), 'DisplayName', 'ode: x2(x1)');
                    title("x2(x1) optimal step 0.2 x2(x1) + ode");
                    legend('show');
                    %
                    figure(5); 
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'step(0.2): x1(t)');
                    title("x1(t) and x2(t) step 0.4 and 0.2 (final and smaller)");
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'step(0.2): x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    title("x2(x1) step 0.4 and 0.2 (final and smaller)")
                    plot(x(1, :), x(2, :), 'DisplayName', 'step(0.2): x2(x1)');
                    legend('show');
                    hold on;
                    %
                end
                figure(1);
                subplot(2,2,1);
                plot(t, x(1,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,2);
                plot(t, x(2,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,3);
                plot(t, x(1,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,4);
                plot(t, x(2,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                figure(2);
                hold on;
                grid on;
                plot(x(1, :), x(2, :), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                legend('show');
                obj.h = obj.h/2;
            end
        end
        function TASK2aP5EC5E(obj)
            obj.h = 0.8;
            figure(1);
            title("Runge-Kutta 4th order, constant step");
            
            subplot(2,2,1);
            hold on;
            grid on;
            title("x1(t)");
            xlim([0, 20]);
            
            subplot(2,2,2);
            hold on;
            grid on;
            title("x2(t)");
            xlim([0, 20]);
            
            subplot(2,2,3);
            hold on;
            grid on;
            title("x1(t) zoom");
            xlim([12.8, 13.4]);
            ylim([0.69, 0.72]);
            
            subplot(2,2,4);
            hold on;
            grid on;
            title("x2(t) zoom");
            xlim([11.3, 11.8]);
            ylim([0.66, 0.72]);

            legend('show');
            for i = 1:1:10
                [t, x] = P5EC5E(obj, obj.h);
                if obj.h == 0.4
                    figure(5); 
                    title("x1(t) and x2(t) larger step size 0.4");
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'step(0.4): x1(t)');
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'step(0.4): x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    title("x2(x1) step 0.4 and 0.2 (final and bigger)")
                    plot(x(1, :), x(2, :), 'DisplayName', 'step(0.4): x2(x1)');
                    legend('show');
                    hold on;
                end
                if obj.h == 0.2
                    [t_ode, x_ode] = ode45(@method, obj.interval, obj.x_val);
                    figure(4);
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'x1(t)');
                    title("x1(t) and x2(t) optimal step 0.2 graph + ode");
                    hold on;
                    plot(t_ode, x_ode(:, 1), 'DisplayName', 'ode: x1(t)');
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'x2(t)');
                    hold on;
                    plot(t_ode, x_ode(:, 1), 'DisplayName', 'ode: x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    plot(x(1, :), x(2, :), 'DisplayName', 'x2(x1)');
                    hold on;
                    plot(x_ode(:, 1), x_ode(:, 2), 'DisplayName', 'ode: x2(x1)');
                    title("x2(x1) optimal step 0.2 x2(x1) + ode");
                    legend('show');
                    %
                    figure(5); 
                    subplot(1, 2, 1);
                    plot(t, x(1,:), 'DisplayName', 'step(0.2): x1(t)');
                    title("x1(t) and x2(t) step 0.4 and 0.2 (final and bigger)");
                    hold on;
                    plot(t, x(2,:), 'DisplayName', 'step(0.2): x2(t)');
                    legend('show');
                    subplot(1, 2, 2);
                    title("x2(x1) step 0.4 and 0.2 (final and bigger)")
                    plot(x(1, :), x(2, :), 'DisplayName', 'step(0.2): x2(x1)');
                    legend('show');
                    hold on;
                    %
                end
                figure(1);
                subplot(2,2,1);
                plot(t, x(1,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,2);
                plot(t, x(2,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,3);
                plot(t, x(1,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                subplot(2,2,4);
                plot(t, x(2,:), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                figure(2);
                hold on;
                grid on;
                plot(x(1, :), x(2, :), 'DisplayName', sprintf("Current step = %0.5f", obj.h));
                legend('show');
                obj.h = obj.h/2;
            end
        end
        function TASK2bRK4(obj)
            obj.h = 4;
            [x, t, error, h] = Runge_Kutta_Variable(obj, 1e-10, 1e-10, 1e-10);
            [t_ode, x_ode] = ode45(@method, obj.interval, obj.x_val);
            hold on;
            figure(1)
            subplot(1,2,1);
            %plot of a generater x2(x1)
            plot(x(1, :), x(2, :));
            hold on;
            plot(x_ode(:, 1), x_ode(:, 2), 'DisplayName', 'ode45 method: x2(x1)');
            legend('show');
            grid on;
            box off;
            title("x2(x1)");
            subplot(1,2,2);
            hold on;
            %plot of a generater x1(t)
            plot(t, x(1, :), 'DisplayName', 'x1(t)');
            hold on;
            plot(t_ode, x_ode(:, 1), 'DisplayName', 'ode45 method: x1(t)');
            %plot of a generater x2(t)
            plot(t, x(2, :), 'DisplayName', 'x2(t)');
            hold on;
            plot(t_ode, x_ode(:, 2), 'DisplayName', 'ode45 method: x2(t)');
            grid on;
            box off;
            title("x1(t) and x2(t)");
            legend('show');
            figure(2)
            subplot(1,2,1);
            %plot of a step
            plot(t, h);
            grid on;
            box off;
            title("step(t)");
            subplot(1,2,2);
            hold on;
            %plot of an error
            plot(t, error);
            grid on;
            box off;
            title("error(t)");
        end
        function ODE(obj, fig)
            figure(fig);
            title('ODE graph');
            [t, x] = ode45(@method, obj.interval, obj.x_val);
            subplot(1,2,1);
            %ODE x2(x1)
            plot(x(:, 1), x(:, 2), 'DisplayName', 'x2(x1)');
            grid on;
            box off;
            title("x2(x1) ODE graph");
            subplot(1, 2, 2);
            %ODE x1(t)
            hold on;
            grid on;
            legend('show');
            plot(t, x(:, 1), 'DisplayName', 'x1(t)');
            title("x1(t) and x2(t) ODE graph");
            %ODE x2(t)
            hold on;
            grid on;
            plot(t, x(:, 2), 'DisplayName', 'x2(t)');
        end
        function [x, t, errors, h] = Runge_Kutta_Variable(obj, error_rel, error_abs, hmin)
            [t, x] = init(obj, obj.h);
            T = 5;
            step = 0.9;
            h(1) = obj.h;
            errors(:, 1) = [0, 0];
            n = 1;

            while t(n) + h(n) < obj.interval(2)
               %this part is calculating the solution
               %basically the same thing as before
               
               k = rk4_in(x, t, h(n), n);
               full_step = x(:, n) + h(n)/6*(k(:, 1)+2*k(:, 2)+2*k(:, 3)+k(:, 4));
               
               x(:, n+1) = full_step;
               %in this part i am calcualting the half-step
               %see page 9 of the report.
               h_2 = h(n)/2;
               half_step = x(:, n);
               for i = 1:1:2
                   k(:, 1) = method(t(n), half_step);
                   k(:, 2) = method(t(n), half_step+h_2*k(:, 1)/2);
                   k(:, 3) = method(t(n), half_step+h_2*k(:, 2)/2);
                   k(:, 4) = method(t(n), half_step+h_2*k(:, 3));
                   
                   half_step = half_step + h_2/6*(k(:, 1)+2*k(:, 2)+2*k(:, 3)+k(:, 4));
               end
               %estimates the error
               %see page 9 of the report
               err = (half_step-full_step)/(2^4-1);
               errors(:, n+1) = err;
               %this is the correction of the step
               %see page 9
               ei = abs(half_step)*error_rel+error_abs;
               %this is on the flow char page 10
               alfa = min((ei./abs(err)).^(1/5));
               h_p = step*alfa*h(n);
               %this part of the code is a representation of the flowchart
               %this is on the flow char page 10
               if step*alfa >= 1
                   if t(n)+h(n) == obj.interval(2)
                       return
                   end
                   t(n+1) = t(n)+h(n);
                   h(n+1) = min([h_p,T*h(n),obj.interval(2)-t(n)]);
                   n = n+1;
               else
                   if h_p < hmin
                       error("Not possible with this hmin");
                   else
                       h(n) = h_p;
                   end
               end
            end
        end
        %The classic RK4 algorithm.
        function [x, t] = Runge_Kutta_Constant(obj, h)
            [t, x] = init(obj, h);
            index = 1;
            %the list of steps can be found on page 8 of the report
            for i = obj.interval(1)+h:h:obj.interval(2)
                k = rk4_in(x, t, h, index);
               
                index = index + 1;
                t(index) = i;
                %final step
                x(:, index) = x(:, index-1) + h/6*(k(:, 1)+2*k(:, 2)+2*k(:, 3)+k(:, 4));
            end
        end
        function [t, x] = P5EC5E(obj, h)
            [t, x] = init(obj, h);
            index = 5;
            %basic parameters
            %for the more information see page 11 of the report
            B_explicit = [1901, -2774, 2616, -1274, 251]/720;
            B_implicit = [475, 1427, -798, 482, -173, 27]/1440;
            
            h_curr = h;
            %7.2.3. Stability and convergence
            %Let up run the RK4 for the first 5 values of x to reduce the error.
            %First value is specified in the init and 1-k are here.
            %If I understand this correctly this approach reduces the error.
            for i = 2:5
                t(i) = h_curr;
                h_curr = h_curr + h;
                k = rk4_in(x, t, h, i-1);
                x(:, i) = x(:, i-1) + h/6*(k(:, 1)+2*k(:, 2)+2*k(:, 3)+k(:, 4));
            end
            
            for i = h_curr:h:obj.interval(2)
                index = index+1;
                t(index) = i;
                %Prediction step
                current_sum = 0;
                for j = 1:1:5
                    if index - j < 1
                        break
                    end
                    current_sum = current_sum + B_explicit(j)*method(t(index-j), x(:, index-j));
                end
                xstep = x(:, index-1) + h*current_sum;
                %Correction
                current_sum = 0;
                for j = 1:1:5
                    if index - j < 1
                        break
                    end
                    current_sum = current_sum + B_implicit(j+1)*method(t(index-j), x(:, index-j));
                end
                x(:, index) = x(:, index-1) + h*current_sum + h*B_implicit(1)*method(t(index-1), xstep);
            end
        end
        function [t, x] = init(obj, h)
            steps = floor(abs(obj.interval(2) - obj.interval(1))/abs(h));
            t = zeros(1, steps+1);
            t(1) = obj.interval(1);
            x(:, 1) = obj.x_val;
        end
    end
end
function k = rk4_in(x, t, h, index)
t(index)
    %k1
    k(:, 1) = method(t(index), x(:, index));
    %k2
    k(:, 2) = method(t(index), x(:, index)+h*k(:, 1)/2);
    %k3
    k(:, 3) = method(t(index), x(:, index)+h*k(:, 2)/2);
    %k4
    k(:, 4) = method(t(index), x(:, index)+h*k(:, 3));
end
function [out] = method(t, x)
    out = [x(2)+x(1)*(0.5-x(1)^2-x(2)^2); -x(1)+x(2)*(0.5-x(1)^2-x(2)^2)];
end