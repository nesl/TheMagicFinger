classdef GpsReceiver < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        measurement_buffer = []
        measurement_ids = [];
        LIGHTSPEED = 299792458; % m/s
    end
    
    methods
        % CONSTRUCTOR
        function obj = GpsReceiver()
        end
        
        % add a gps measurement if it's unique
        function addMeasurement(obj, meas)
            %if isempty( find(obj.measurement_ids == meas.id ) )
                obj.measurement_buffer = [obj.measurement_buffer; meas];
                obj.measurement_ids = [obj.measurement_ids; meas.id];
            %end
        end
        
        % get estimated position from measurements
        function p = estimatePosition(obj, x_init)
            % construct data vector
            %{
            X = [];
            for i=1:length(obj.measurement_buffer)
                meas = obj.measurement_buffer(i);
                X = [X; meas.vectorize()];
            end
            %}
            
            % run nonlinear least squares
            fcn = @(x) obj.costFunction(x);
            %options = optimoptions('lsqnonlin','TolFun',0.1e-9);
            options = optimoptions('lsqnonlin','Display','none', 'Algorithm', 'trust-region-reflective', 'TolX', 1e-9, 'TolFun', 1e-9);
            [p,resnorm] = lsqnonlin(fcn, x_init, [], [], options);
            fprintf('resnorm = %.3f\n', resnorm);
        end
        
        % cost function to reduce for position estimate 
        % NEED TO VECTORIZE SO IT IS SUM OF SQUARES !!! OTHERWISE WILL ONLY
        % CONSTRAIN ONE VARIABLE, E.G. THE OFFSET !!!
        function C = costFunction(obj, x)
            % empty cost vector to start
            C = zeros(1,length(obj.measurement_buffer));
            
            % starting local time of measurements in this group
            tbuffer_start = obj.measurement_buffer(1).time_rx;
            
            % extract state from vector
            offset = x(1);
            bias = x(2);
            xyz = x(3:5);
                        
            % add cost for each measurement
            for i=1:length(obj.measurement_buffer)
                % extract message
                m = obj.measurement_buffer(i);
                % time since first message in the buffer
                dt = m.time_rx - tbuffer_start;
                % range & prop. delay to this hypothetical position
                range = norm( m.xyz - xyz );
                Tp = range/obj.LIGHTSPEED;
                % convert local rx time to global time
                global_rx = m.time_rx - offset - dt*bias;
                % compare predicted and measured prop. delay
                global_Tp = global_rx - m.globalTxTime();
                mcost = (Tp - global_Tp)*1e9; % ns for precision
                % append to cost vector
                C(i) = mcost;
            end
        end
        
        % clear measurement buffer
        function clearBuffer(obj)
            obj.measurement_buffer = [];
            obj.measurement_ids = [];
        end
        
        % how many unique messages are there?
        function n = numUniqueMessages(obj)
            n =  length( unique(obj.measurement_ids) );
        end
        
        % how many total messages are there?
        function n = numTotalMessages(obj)
            n = length(obj.measurement_buffer);
        end
    end
    
end

