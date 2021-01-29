classdef TASK1
    properties
        x; %x initial points
        y; %y initial points
    end
    
    methods
        function obj = TASK1()
            obj.x = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5];
            obj.y = [-7.2606, -5.6804, -3.7699, -3.4666, -1.8764, -2.1971, -2.8303, -6.3483, -13.3280, -23.4417, -36.2458];
        end
        function showInitialPoints(obj)
            plot(obj.x, obj.y, 'ob');
            grid on;
            title("Initial points");
            xlim([-5, 5]);
        end
        %Executes an LS_App_WQR algorithm for degree.
        function ExecuteForDegrees(obj, degree)
            figure(1)
            plot(obj.x, obj.y, 'ob', 'DisplayName', 'Initial points');
            grid on;
            title("Least Square Approximation");
            hold on;
            legend('show', 'Location', 'southwest');
            PLOTx = linspace(-5, 5);
            error = zeros(1, degree+1);
            polynomial = zeros(degree+1);
            for curr = 0:1:degree
                %Approximiation the individual degree.
                [ind, cond] = LS_App_WQR(obj, curr);
                fprintf("Condition number of Gram's: %d, degree: %d \n", cond, curr);
                %Adding to the plot.
                PLOTy = polyval(ind, PLOTx);
                plot(PLOTx, PLOTy, 'DisplayName', sprintf("Degree = %d", curr));
                %--------------------
                %Errors calculation.
                solutions = polyval(ind, obj.x);
                %Error and the approximation.
                error(curr+1) = norm(obj.y - solutions);
                
            end
            hold off;
            figure(2);
            title("Least Square Approximation My Solution");
            grid on;
            hold on; 
            plot(obj.x, obj.y, 'ob', 'DisplayName', 'Initial points');
            legend('show', 'Location', 'southwest');
            plot(PLOTx, PLOTy, 'DisplayName', sprintf("Degree = %d", curr));
            hold off;
            %calculates the results in a form of actual f(x) function.
            [n, m] = size(polynomial);
            for i = 1:1:n
                fprintf("f(x) = ");
                for j = 1:1:m
                    if polynomial(i, j) ~= 0
                        if(j == degree + 1)
                            fprintf("%f\n", polynomial(i, j));   
                        else
                            if polynomial(i, j) < 0
                                if m-j ~= 1
                                    fprintf("%fx^%d ", polynomial(i, j), m-j);
                                else
                                    fprintf("%fx ", polynomial(i, j));
                                end
                            else
                                if m-j ~= 1
                                    fprintf("+%fx^%d ", polynomial(i, j), m-j);
                                else
                                    fprintf("+%fx ", polynomial(i, j));
                                end
                            end
                        end
                    end
                end
            end
            %norm
            for i = 1:1:length(error)
                fprintf("Norm for power: %d -> %.10f\n", i - 1, error(i));
            end
        end
        function [solution, condN] = LS_App_WQR(obj, degree)
            %A method finds a polynomial of a given degree.
            %Which fits the given set of points obj.x and obj.y.
            %A least squares aprroximation with QR factorization is used.
            SIZE = size(obj.y);
            if SIZE(1) == 1
                obj.y = obj.y';
            end
            out = zeros(length(obj.x), degree+1);
            for i = 1:1:length(obj.x)
                for j = 0:1:degree
                    out(i, j+1) = obj.x(i)^j;
                end
            end
            condN = cond(out' * out);
            %QR factorization of the out matrix.
            [Q, R] = QRFact(out);
            solution = fliplr((R\(Q'*obj.y))');
        end
        
    end
end
function [Q,R] = QRFact(in)
    %QR Factorization of a matrix A
    %outputs the Q and R
    [m, n] = size(in);
    Q = zeros(m, n);
    R = zeros(n, n);
    d = zeros(1, n);
    %Factorization of an input.
    for i = 1:1:n
        Q(:, i) = in(:, i);
        R(i, i) = 1;
        d(i) = Q(:, i)' * Q(:, i);
        for j = i+1:1:n
            R(i, j) = (Q(: ,i)' * in(:, j)) / d(i);
            in(:, j) = in(:, j) - R(i, j) * Q(:, i);
        end
    end
    %Normalization of an input.
    for i = 1:1:n
        dd = norm(Q(:, i));
        Q(:, i) = Q(:, i) / dd;
        R(i, i:n) = R(i, i:n) * dd;
    end
end
