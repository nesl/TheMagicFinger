classdef GpsMessage < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id = 0;
        xyz = [0 0 0];
        time_tx = 0;
        time_rx = 0;
        offset = 0;
        bias
    end
    
    methods
        function obj = GpsMessage(anc_id, xyz, time_tx, time_rx, anc_ofst, anc_bias)
            obj.id = anc_id;
            obj.xyz = xyz;
            obj.time_tx = time_tx;
            obj.time_rx = time_rx;
            obj.offset = anc_ofst;
            obj.bias = anc_bias;
        end
        
        % transform a local, rx time to global frame
        function g = globalTxTime(obj)
            g = obj.time_tx - obj.offset;
        end
        
        % get the anchor's location
        function xyz = getAnchorPosition(obj)
            xyz = self.xyz;
        end
        
        % vectorize this message
        function v = vectorize(obj)
            v = [ 
                obj.xyz;
                obj.transformTime( obj.time_tx )
                ];
        end
    end
    
end

