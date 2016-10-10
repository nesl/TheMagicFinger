classdef TdoaServer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        measurement_buffer = []
        measurement_seqs = [];
        LIGHTSPEED = 299792458; % m/s
    end
    
    methods
        % CONSTRUCTOR
        function obj = TdoaServer()
        end
        
        % add a gps measurement if it's unique
        function addMeasurement(obj, meas)
            obj.measurement_buffer = [obj.measurement_buffer; meas];
            obj.measurement_seqs = [obj.measurement_seqs; meas.msg_seq];
        end
        
        % get measurements with a specific sequence number
        function measList = getMeasurementsWithSeq(obj, seq)
            measList = [];
            for i=1:length(obj.measurement_buffer)
                if obj.measurement_buffer(i).msg_seq == seq
                    measList = [measList obj.measurement_buffer(i)];
                end
            end 
        end
        
        % get estimated position from measurements
        function p = estimatePosition(obj, seq, x_init)
            % collect sequence measurements
            seqMeasurements = obj.getMeasurementsWithSeq(seq);
            
            % convert into global RX times
            trx = [];
            ids = [];
            for i=1:length(seqMeasurements)
                trx = [trx; (seqMeasurements(i).time_rx - seqMeasurements(i).anc_ofst)];
                ids = [ids; seqMeasurements(i).anc_id];
            end
            
            fprintf('--------------------\n');
            for i=1:length(seqMeasurements)
                fprintf('[%d]: TRX = %.3f\n', ids(i), (trx(i) - min(trx))*1e9);
            end
            
            % nonlinear least squares
            fcn = @(x) obj.costFunction(seqMeasurements, x);
            options = optimoptions('lsqnonlin','Display','none', 'Algorithm', 'trust-region-reflective', 'TolX', 1e-9, 'TolFun', 1e-9);
            [p,resnorm] = lsqnonlin(fcn, x_init, [], [], options);
        end
        
        % cost function to reduce for position estimate 
        % NEED TO VECTORIZE SO IT IS SUM OF SQUARES !!! OTHERWISE WILL ONLY
        % CONSTRAIN ONE VARIABLE, E.G. THE OFFSET !!!
        function C = costFunction(obj, measList, xyz)
            % empty cost vector to start
            C = zeros(1,length(measList));
            
            % add cost for each measurement
            for i=1:length(measList)
                for j=i:length(measList)
                    if j==i
                        continue;
                    end
                    % pair of anchors
                    m1 = measList(i);
                    m2 = measList(j);
                    
                    % receive times in global frame
                    gRx1 = m1.time_rx - m1.anc_ofst;
                    gRx2 = m2.time_rx - m2.anc_ofst;
                    
                    % expected propagation delays
                    range1 = norm( m1.anc_xyz - xyz );
                    range2 = norm( m2.anc_xyz - xyz );
                    
                    % range difference
                    rangeDiff = range1 - range2;
                    
                    % difference in measured propagation distance
                    measDiff = (gRx1 - gRx2)*obj.LIGHTSPEED;

                    % append to cost vector
                    C(i) = (rangeDiff - measDiff);
                end
            end
        end
        
        % clear measurement buffer
        function clearBuffer(obj)
            obj.measurement_buffer = [];
            obj.measurement_seqs = [];
        end
        
        % how many unique messages are there?
        function [seq,num] = getSeqCounts(obj)
            if isempty(obj.measurement_buffer)
                seq = 0;
                num = 0;
                return;
            end
            [m,f] = mode(obj.measurement_seqs);
            seq = m(1);
            num = f(1);
        end
        
    end
    
end

