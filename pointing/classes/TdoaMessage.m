classdef TdoaMessage < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        target_id = 0;
        anc_id = 0;
        msg_seq = 0;
        anc_xyz = [0 0 0];
        anc_ofst = 0;
        time_tx = 0;
        time_rx = 0;
    end
    
    methods
        function obj = TdoaMessage(target_id, anc_id, msg_seq,  anc_xyz, anc_ofst, time_tx, time_rx)
            obj.target_id = target_id;
            obj.anc_id = anc_id;
            obj.msg_seq = msg_seq;
            obj.anc_xyz = anc_xyz;
            obj.anc_ofst = anc_ofst;
            obj.time_tx = time_tx;
            obj.time_rx = time_rx;
        end
        
        % transform a local, rx time to global frame
        function g = globalRxTime(obj)
            g = obj.time_rx - obj.anc_ofst;
        end
        
        % get the anchor's location
        function xyz = getAnchorPosition(obj)
            xyz = self.xyz;
        end
    end
    
end

