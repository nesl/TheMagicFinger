classdef IotManager < handle
    %NETWORKMANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        LIGHTSPEED = 299792458; % m/s
        PLTSTYLES = {'-o', '-s', '-^', '-x', '-*', '--o', '--s', '--^', '--x', '--*'};
        PLTCOLORS = [
            1 0 0
            0 0 1
            0 0 0
            0.8 0.8 0
            1 0 1
            0 0.7 0.7
            0 1 0
            0.5 0.5 0.5
            0.25 0.25 1
            1 0.25 0.25
            0.25 1 0.25
            ];
        PASS = 1;
        BLOCK = 0;
        MSGTYPE1 = 1;
        MSGTYPE2 = 2;
        MSGTYPE3 = 3;
        MEASTYPE_d = 1;
        MEASTYPE_r = 2;
        MEASTYPE_b = 3;
        MEASTYPE_B = 4;
        MEASTYPE_R = 5;
        MAXMEASUREMENTS = 1e6;
        MEASBLOCKBKOFF = 0.010;
    end
    
    properties
        % node variables
        nodes = {};
        nodeinfo = {};
        numnodes = 0;

        % message / data variables
        dataparser = [];
        antennacorrections = [];
        twrcorrections = [];
        owrcorrections = [];
        msgidx = 1;
        measurement_history = {};
        meascnt = 0;
        
        % range bias state
        rangebiases = [];
        var_rangebias = 0.10;  % m/s
        vari_rangebias = 1.0^2; % m
        
        % message filtering (pass=1,block=0)
        filter_type1 = 0;
        filter_type2 = 0;
        filter_type3 = 0;
        filter_meas_d = 0;
        filter_meas_r = 0;
        filter_meas_b = 0;
        filter_meas_B = 0;
        filter_meas_R = 0;
        
        % min period between message types (sec)
        period_type1 = 0;
        period_type2 = 0;
        period_type3 = 0;
        
        % node blacklist
        blacklist = [];
        
        % message queues (src, dest, type)
        measQ = [];
        
        % reference node id
        refNode = [];
        
        % enable SLATS?
        enable_slats = false;
    end
    
    methods
        % =============================================
        %                  CONSTRUCTOR
        % =============================================
        function obj = IotManager(configfile, logfolder, varargin)
            % create dataparser
            obj.dataparser = DataParserROS(configfile, logfolder);
            obj.nodeinfo = obj.dataparser.getNodeInfo();
            % create node objects
            for i=1:length(obj.nodeinfo)
                node_names = obj.nodeinfo{1};
                node_ids = obj.nodeinfo{2};
                node_dynamics = obj.nodeinfo{6};
                for i=1:length(node_ids)
                    name = node_names{i};
                    id = node_ids(i);
                    xyz = [obj.nodeinfo{3}(i); obj.nodeinfo{4}(i); obj.nodeinfo{5}(i)];
                    obj.nodes{i} = Node(id, name);
                    obj.nodes{i}.setTruePosition(xyz);
                    % is this a mobile node?
                    if strcmp(node_dynamics{i}, 'mobile');
                        obj.nodes{i}.setToMobile();
                    end
                end
            end
            obj.numnodes = length(obj.nodeinfo{1});
            
            % load rangefile corrections if they exist % TWR, OWR
            if nargin >= 3
                ant_file = varargin{1};
                obj.antennacorrections = csvread(ant_file);
            end
            if nargin >= 4
                owr_file = varargin{2};
                obj.owrcorrections = csvread(owr_file);
            end
            if nargin >= 5
                twr_file = varargin{3};
                obj.twrcorrections = csvread(twr_file);
            end
            
            % pre-allocate measurement memory
            obj.measurement_history = cell(obj.MAXMEASUREMENTS,1);
            
            % pre-allocate range bias memory
            obj.rangebiases = zeros( obj.getNumNodes(), obj.getNumNodes() );
        end
        
        % =============================================
        %                MESSAGE FILTERING
        % =============================================
        
        % queue a specific message type
        function queueMessage(obj, srcId, dstId, type)
            % add if Q empty
            if isempty(obj.measQ)
                obj.measQ = [obj.measQ; [srcId dstId type]];
                return;
            end
            % otherwise only add if not already queued
            if isempty( find( obj.measQ(:,1) == srcId & obj.measQ(:,2) == dstId & obj.measQ(:,3) == type ) )
                obj.measQ = [obj.measQ; [srcId dstId type]];
            end
        end
        
        % enable / disable specific message types
        function enableMessageType(obj, type, enable)
            switch type
                case obj.MSGTYPE1
                    obj.filter_type1 = enable;
                case obj.MSGTYPE2
                    obj.filter_type2 = enable;
                case obj.MSGTYPE3
                    obj.filter_type3 = enable;
                otherwise
                    error('Message type not recognized');
            end
        end
        
        % downsample specific message types to a given period lower bound
        function setMessageMinInterval(obj, type, interval)
            switch type
                case obj.MSGTYPE1
                    obj.period_type1 = interval;
                case obj.MSGTYPE2
                    obj.period_type2 = interval;
                case obj.MSGTYPE3
                    obj.period_type3 = interval;
                otherwise
                    error('Message type not recognized');
            end
        end
        
        % enable / disable specific message types
        function enableMeasurementType(obj, type, enable)
            switch type
                case obj.MEASTYPE_d
                    obj.filter_meas_d = enable;
                case obj.MEASTYPE_r
                    obj.filter_meas_r = enable;
                case obj.MEASTYPE_b
                    obj.filter_meas_b = enable;
                case obj.MEASTYPE_B
                    obj.filter_meas_B = enable;
                case obj.MEASTYPE_R
                    obj.filter_meas_R = enable;
                otherwise
                    error('Message type not recognized');
            end
        end
        
        % determine whether or not this measurement should be received
        function valid = filterMeasurement(obj, meas)
            srcId = meas.getSourceId();
            dstId = meas.getDestId();
            type = meas.getType();
            tlastSent = obj.nodes{ obj.getNodeIdx(srcId) }.getLastMeasTime( type );
            
            % ignore measurement if src or dst are blacklisted
            if ~isempty( find(obj.blacklist == srcId, 1) ) || ~isempty( find(obj.blacklist == dstId, 1) )
                valid = obj.BLOCK;
                return;
            end
            
            % set allowed measurement types
            if obj.filter_meas_d
                meas.allowMeasType_d();
            end
            if obj.filter_meas_r
                meas.allowMeasType_r();
            end
            if obj.filter_meas_R
                meas.allowMeasType_R();
            end
            
            % next, pass immediately if measurement type has been requested
            if ~isempty( obj.measQ )
                qIdx = find(obj.measQ(:,1) ==  srcId & obj.measQ(:,2) == dstId);
                if ~isempty(qIdx)
                    % force this message to be the right type
                    meas.type = obj.measQ(qIdx,3);
                    meas.queued = 1;
                    % a message of this type has been requested, pass it and remove
                    % from the message queue
                    
                    % block if the range or time is obviously incorrect
                    if meas.r_ij < -1 || meas.r_ij > 20 || meas.R_ij < -1 || meas.R_ij > 20
                        valid = obj.BLOCK;
                        obj.measQ(qIdx,:) = [];
                        return;
                    end
                    valid = obj.PASS;
                    obj.measQ(qIdx,:) = [];
                    return;
                end
            end
            
            % block if no measurements supported
            if ~(obj.filter_meas_d || obj.filter_meas_r || obj.filter_meas_R)
                valid = obj.BLOCK;
                return;
            end
                        
            % block if the range or time is obviously incorrect
            if meas.r_ij < -1 || meas.r_ij > 20 || meas.R_ij < -1 || meas.R_ij > 20
                valid = obj.BLOCK;
                return;
            end
              
            
            % pass as a type 3 message if they're allowed and it's been
            % long enough since the last type 3 message
            if obj.filter_type3 == obj.PASS
                % downgrade this to a type 2 message
                meas.setType(obj.MSGTYPE3);
                tlastSent = obj.nodes{ obj.getNodeIdx(srcId) }.getLastMeasTime( obj.MSGTYPE3 );
                dt = meas.getTime() - tlastSent;
                if dt >= obj.period_type3 || dt < obj.MEASBLOCKBKOFF
                    valid = obj.PASS;
                    obj.nodes{ obj.getNodeIdx(srcId) }.setLastMeasTime( obj.MSGTYPE3, meas.getTime() );
                    return;
                end
            end
            % otherwise make this a type 2 message and pass it if they're
            % allowed and if it's been long enough since the last type 2.
            if obj.filter_type2 == obj.PASS
                % downgrade this to a type 2 message
                meas.setType(obj.MSGTYPE2);
                tlastSent = obj.nodes{ obj.getNodeIdx(srcId) }.getLastMeasTime( obj.MSGTYPE2 );
                dt = meas.getTime() - tlastSent;
                if dt >= obj.period_type2 || dt < obj.MEASBLOCKBKOFF
                    valid = obj.PASS;
                    obj.nodes{ obj.getNodeIdx(srcId) }.setLastMeasTime( obj.MSGTYPE2, meas.getTime() );
                    return;
                end
            end
            % finally, we can at the very least make this a type 1 message
            % if the other 2 messages are blocked for this time period.
            % Disallow this too, though, if not enough time has elapsed
            % since the last type 1.
            if obj.filter_type1 == obj.PASS
                % downgrade this to a type 1 message
                meas.setType(obj.MSGTYPE1);
                tlastSent = obj.nodes{ obj.getNodeIdx(srcId) }.getLastMeasTime( obj.MSGTYPE1 );
                dt = meas.getTime() - tlastSent;
                if dt >= obj.period_type1 || dt < obj.MEASBLOCKBKOFF
                    valid = obj.PASS;
                    obj.nodes{ obj.getNodeIdx(srcId) }.setLastMeasTime( obj.MSGTYPE1, meas.getTime() );
                    return;
                end
            end
            
            % otherwise, block it
            valid = obj.BLOCK;
        end
        
        % =============================================
        %                  PLOTTING
        % =============================================
        
        % plot styles and colors
        function style = getPlotStyle(obj, nodeId)
            nidx = obj.dataparser.getNodeIdx(nodeId);
            style = obj.PLTSTYLES{nidx};
        end
        function c = getPlotColor(obj, nodeId)
            nidx = obj.dataparser.getNodeIdx(nodeId);
            c = obj.PLTCOLORS(nidx,:);
        end
        
        % get number of nodes
        function n = getNumNodes(obj)
            n = obj.numnodes;
        end
        
        % get number of static nodes
        function n = getNumStaticNodes(obj)
            n = 0;
            for i=1:obj.getNumNodes()
                dynamics = obj.nodeinfo{6}{i};
                if strcmp(dynamics, 'static')
                    n = n + 1;
                end
            end
        end
        
        % get number of mobile nodes
        function n = getNumMobileNodes(obj)
            n = obj.getNumNodes() - obj.getNumStaticNodes();
        end
        
        % =============================================
        %                  NODE INFO
        % =============================================
        
        % get node ascii names 
        function n = getNodeNames(obj)
            n = obj.nodeinfo{1};
        end
        
        % get numeric node is: 0, 1, ...
        function id = getNodeIds(obj)
            id = obj.nodeinfo{2};
        end
        
        % get the index of a specific node
        function idx = getNodeIdx(obj, nodeid)
            idx = obj.dataparser.getNodeIdx( nodeid );
        end
        
        % blacklist a certain node (don't use it)
        function blacklistNode(obj, nodeid)
            idx = obj.getNodeIdx( nodeid );
            nodeids = obj.getNodeIds();
            id = nodeids(idx);
            obj.blacklist = [obj.blacklist; id];
        end
        
        % =============================================
        %               STATE ACCESSORS
        % =============================================
        
        % get network state
        function s = getState(obj)
            s = [];
            for i=1:length(obj.nodes)
                s = [s; obj.nodes{i}.getState()];
            end
        end
        
        % get the vectorized range biases
        function sb = getRangeBiasVec(obj)
            sb = [];
            for i=1:obj.getNumNodes()
                for j=(i+1):obj.getNumNodes()
                    sb = [sb; obj.rangebiases(i,j)];
                end
            end
        end
        
        % range bias vector to matrix
        function SB = rangeBiasVec2Mat(obj, sb)
            SB = [];
            idx = 0;
            for i=1:obj.getNumNodes()
                for j=(i+1):obj.getNumNodes()
                    idx = idx + 1;
                    SB(i,j) = sb(idx);
                    SB(j,i) = sb(idx);
                end
            end
        end
        
        % set the vectorized range biases
        function setRangeBiasState(obj, sb)
            obj.rangebiases = obj.rangeBiasVec2Mat(sb);
        end
        
        % set network state
        function setState(obj, s)
            idx = 1;
            for i=1:length(obj.nodes)
                stateSize = length(obj.nodes{i}.getState());
                si = s(idx:(idx+stateSize-1));
                obj.nodes{i}.setState(si);
                idx = idx + stateSize;
            end

        end
        
        % get estimated node positions for only static nodes
        function P = getEstimatedStaticPositions(obj)
            P = [];
            for i=1:length(obj.nodes)
                if ~obj.nodes{i}.isMobile()
                    P = [P; obj.nodes{i}.getId() obj.nodes{i}.getStatePosition()'];
                end
            end
        end
        
        % get estimated node positions for all nodes
        function P = getEstimatedPositions(obj)
            P = [];
            for i=1:length(obj.nodes)
                P = [P; obj.nodes{i}.getStatePosition()'];
            end
        end
        
        % get true node positions
        function P = getTrueStaticPositions(obj)
            P = [];
            for i=1:length(obj.nodes)
                if ~obj.nodes{i}.isMobile()
                    P = [P; obj.nodes{i}.getId() obj.nodes{i}.getTruePosition()'];
                end
            end
        end
        
        % get estimated positions transformed to match true as closely as
        % possible (procrustes)
        function P = getTransformedPositions(obj)
            truStatic = obj.getTrueStaticPositions();
            estStatic = obj.getEstimatedStaticPositions();
            % procrustes ignoring first column (node IDs)
            [~,~,Transform] = procrustes(truStatic(:,2:end), estStatic(:,2:end), 'Scaling', false);
            % transform all points
            estAll = obj.getEstimatedPositions();
            ofst = repmat(Transform.c(1,:), obj.getNumNodes(), 1);
            P = Transform.b*estAll*Transform.T + ofst;
        end
        
        % set reference node (for EKF)
        function setReferenceNode(obj, nodeId)
            nidx = obj.dataparser.getNodeIdx(nodeId);
            obj.nodes{nidx}.setAsReference();
            obj.refNode = obj.nodes{nidx}.getId();
        end
        
        % request the next (filtered) measurement
        function m = getNextMeasurement(obj)
            while obj.msgidx < obj.dataparser.getNumMeasurements()
                % get new tentative meas
                raw = obj.dataparser.getMeasurement(obj.msgidx);
                meas = Measurement(obj.MSGTYPE3, raw);
                obj.msgidx = obj.msgidx + 1;
                % check if we should pass this message
                if obj.filterMeasurement(meas) == obj.PASS
                    % ANT corrections
                    if ~isempty(obj.antennacorrections)
                        idx = find(obj.antennacorrections(:,1) == meas.getSourceId() &...
                            obj.owrcorrections(:,2) == meas.getDestId());
                        if ~isempty(idx)
                            meas.d_ij = meas.d_ij + 0*1e-9*obj.antennacorrections(idx, 3);
                        end
                    end
                    % OWR corrections
                    if ~isempty(obj.owrcorrections)
                        idx = find(obj.owrcorrections(:,1) == meas.getSourceId() &...
                            obj.owrcorrections(:,2) == meas.getDestId());
                        if ~isempty(idx)    
                            meas.r_ij = meas.r_ij + obj.owrcorrections(idx, 3);
                        end
                    end
                    % TWR corrections
                    if ~isempty(obj.twrcorrections)
                        idx = find(obj.twrcorrections(:,1) == meas.getSourceId() &...
                            obj.twrcorrections(:,2) == meas.getDestId());
                        if ~isempty(idx)    
                            meas.R_ij = meas.R_ij + obj.twrcorrections(idx, 3);
                        end
                    end
                    
                    m = meas;
                    
                    % increment and check message count
                    if obj.meascnt > obj.MAXMEASUREMENTS
                        m = [];
                        return;
                    end
                    obj.meascnt = obj.meascnt + 1;
                    obj.measurement_history{obj.meascnt} = m;
                    return;
                end
            end
            % no more measurements
            m = [];
        end
        
        function killLastMeasurement(obj)
            obj.meascnt = obj.meascnt - 1;
        end
        
        % reset measurement index
        function resetMeasurements(obj)
            obj.msgidx = 1;
            obj.measurement_history = cell(obj.MAXMEASUREMENTS,1);
            obj.meascnt = 0;
        end
        
        % =============================================
        %                  BOOTSTRAPPING
        % =============================================
        
        % skip some time
        function skipTime(obj, tskip)
            % enable type 3 messages for this
            type3Enabled = obj.filter_type3;
            typeREnabled = obj.filter_meas_R;
            obj.enableMessageType(obj.MSGTYPE3, true);
            obj.enableMeasurementType(obj.MEASTYPE_R, true);
           
            % get the first measurement
            meas1 = obj.getNextMeasurement();
            t1 = meas1.getTime();
            t_now = t1;
            while t_now - t1 < tskip
                meas = obj.getNextMeasurement();
                if isempty(meas)
                    break;
                end
                t_now = meas.getTime();
            end
            
            % reset the measurement history but not the dp index
            obj.measurement_history = cell(obj.MAXMEASUREMENTS,1);
            obj.meascnt = 0;
            
            % set type 3 filtering to whatever it was before bootstrap
            obj.enableMessageType(obj.MSGTYPE3, type3Enabled);
            obj.enableMeasurementType(obj.MEASTYPE_R, typeREnabled);
        end
        
        % bootstrap node positions using ranges averaged over a duration
        function bootstrapNodeStates(obj, T_boot)
            % enable type 3 messages for this
            type3Enabled = obj.filter_type3;
            typeREnabled = obj.filter_meas_R;
            obj.enableMessageType(obj.MSGTYPE3, true);
            obj.enableMeasurementType(obj.MEASTYPE_R, true);
            
            % get range estimates (distance matrix)
            R_sum = zeros( obj.getNumNodes(), obj.getNumNodes() );
            % get offset estimates
            D_sum = zeros( obj.getNumNodes(), 1 );
            % distance measurement counts
            D_cnt = zeros( obj.getNumNodes(), obj.getNumNodes() );
            % timing measurements
            d_ij_arrays = cell( obj.getNumNodes(), 1 );
            d_ij_times  = cell( obj.getNumNodes(), 1 );
            % get the first measurement
            meas1 = obj.getNextMeasurement();
            t1 = meas1.getTime();
            t_now = t1;
            while t_now - t1 < T_boot
                meas = obj.getNextMeasurement();
                if isempty(meas)
                    break;
                end
                t_now = meas.getTime();
                srcIdx = obj.getNodeIdx( meas.getSourceId() );
                dstIdx = obj.getNodeIdx( meas.getDestId() );
                
                R_sum(srcIdx,dstIdx) = R_sum(srcIdx,dstIdx) + meas.R_ij;
                R_sum(dstIdx,srcIdx) = R_sum(dstIdx,srcIdx) + meas.R_ij;
                D_cnt(srcIdx,dstIdx) = D_cnt(srcIdx,dstIdx) + 1;
                D_cnt(dstIdx,srcIdx) = D_cnt(dstIdx,srcIdx) + 1;
                
                % look for messages from the reference node to bootstrap
                % clock offset and bias
                if ~isempty(obj.refNode) && meas.getSourceId() == obj.refNode
                    % append d_ij and time
                    d_ij_arrays{dstIdx} = [d_ij_arrays{dstIdx} meas.d_ij];
                    d_ij_times{dstIdx} = [d_ij_times{dstIdx} meas.getTime()];
                end
            end
            
            % average range estimates
            D = zeros( obj.getNumNodes(), obj.getNumNodes() );
            W = ones( obj.getNumNodes(), obj.getNumNodes() );
            for i=1:obj.getNumNodes()
                for j=1:obj.getNumNodes()
                    % ignore diagonals
                    if i == j
                        continue;
                    end
                    % assign default or average measurements
                    if D_cnt(i,j) == 0
                        % if we have no measurement, assume it's a bit far
                        % away and give it very low weight
                        D(i,j) = 10;
                        W(i,j) = 0.01;
                    else
                        D(i,j) = R_sum(i,j) / D_cnt(i,j);
                    end
                end
            end
                        
            % use MDS to get an estimate of the node positions
            Y = mdscale(D, 3, 'Weights', W);
            
            % assign the estimated positions to each node
            for i=1:obj.getNumNodes()
                obj.nodes{i}.setStatePosition( Y(i,:) );
            end
            
            % perform linear regression over d_ij and time to get estimated
            % clock offset and bias w.r.t. reference node
            for i=1:size(d_ij_arrays,1)
                if length(d_ij_arrays{i}) >= 3
                    % wall time = x
                    x = d_ij_times{i} - d_ij_times{i}(1);
                    % measured d_ij = y
                    y = d_ij_arrays{i};
                    p = polyfit(x, y, 1);
                    est_bias = p(1);
                    est_ofst = polyval(p, x(end));
                    obj.nodes{i}.setStateClockOfst(est_ofst);
                    obj.nodes{i}.setStateClockBias(est_bias*1e9);
                    fprintf('BOOTSTRAP: node %d: offset= %.6f ms, bias = %d ppb\n', obj.nodes{i}.getId(), est_ofst*1e3, est_bias*1e9);
                elseif obj.nodes{i}.getId() ~= obj.refNode
                    fprintf('WARN: No timing bootstrap for %d-->%d (%d measurements)\n', obj.refNode, obj.nodes{i}.getId(), length(d_ij_arrays{i}));
                end
            end
            
            % reset the measurement history but not the dp index
            obj.measurement_history = cell(obj.MAXMEASUREMENTS,1);
            obj.meascnt = 0;
            
            % set type 3 filtering to whatever it was before bootstrap
            obj.enableMessageType(obj.MSGTYPE3, type3Enabled);
            obj.enableMeasurementType(obj.MEASTYPE_R, typeREnabled);
        end
        
        % set all static nodes to their true positions
        function setStaticNodesToTruePosition(obj)
            for i=1:length(obj.nodes)
                if ~obj.nodes{i}.isMobile()
                    obj.nodes{i}.setPositionToTrue();
                end
            end
        end
        
        % =============================================
        %              MEASUREMENT HISTORY
        % =============================================
        
        % get pairwise TWR range measurements
        function r = getMeasurementTWR(obj, srcId, dstId)
            r = [];
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                if meas.getSourceId() == srcId && meas.getDestId() == dstId
                    r = [r; meas.R_ij];
                end
            end
        end
        
        % get pairwise OWR range measurements
        function r = getMeasurementOWR(obj, srcId, dstId)
            r = [];
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                if meas.getSourceId() == srcId && meas.getDestId() == dstId
                    r = [r; meas.r_ij];
                end
            end
        end
        
        % get pairwise offset measurements
        function [o,t] = getMeasurementOffsets(obj, srcId, dstId)
            o = [];
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                if meas.getSourceId() == srcId && meas.getDestId() == dstId
                    o = [o; meas.d_ij];
                end
            end
        end
        
        % get all measurements between two devices
        function m = getMeasurements(obj, srcId, dstId)
            m = [];
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                if meas.getSourceId() == srcId && meas.getDestId() == dstId
                    m = [m; meas];
                end
            end
        end
        
        % get pairwise measurement times
        function t = getMeasurementTimes(obj, srcId, dstId)
            t = [];
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                if meas.getSourceId() == srcId && meas.getDestId() == dstId
                    t = [t; meas.getTime()];
                end
            end
        end
        
        % get all measurement times
        function t = getAllMeasurementTimes(obj)
            t = zeros(1,obj.meascnt);
            for i=1:obj.meascnt
                meas = obj.measurement_history{i};
                t(i) = meas.getTime();
            end
        end
        
        % get allan deviation between two nodes
        function [ad,tau] = getAllanDev(obj, srcId, dstId, tau)
            % get non-zero aligned data
            data = obj.dataparser.aligned_logs;
            data = data(data(:,2) == srcId & data(:,3) == dstId, :);
            % get tx and rx
            tx = data(:,5);
            rx = data(:,6);
            % tx time (ignore 1st b/c y is diffed)
            x = tx(2:end);
            % rx dFreq in seconds
            y = (diff(rx) ./ diff(tx) - 1);
            [ad,~,~,tau] = allan(struct('freq',y,'time',x),tau,'',0);
        end
        
        % get all time offsets between two nodes
        function [o,t] = getAllOffsets(obj, srcId, dstId)
            data = obj.dataparser.aligned_logs;
            idxs = find(data(:,2) == srcId & data(:,3) == dstId);
            t = data(idxs,1);
            o = data(idxs,10) - data(idxs,9);
        end
        
        % get all time biases between two nodes
        function [b,t] = getAllBiases(obj, srcId, dstId)
            data = obj.dataparser.aligned_logs;
            idxs = find(data(:,2) == srcId & data(:,3) == dstId);
            t = data(idxs,1);
            o = data(idxs,10) - data(idxs,9);
            b = diff(o)./diff(t);
            t = t(2:end);
        end
        
        % get measurements from a single anchor in a given time period
        function [t,r1] = getRangesInWindow(obj, srcId, dstId, tstart, tstop)
            data = obj.dataparser.timinglogs;
            idxs = find(data(:,2) == srcId & data(:,3) == dstId & data(:,1) >= tstart & data(:,1) <= tstop);
            t = data(idxs,1);
            t0 = data(idxs,5);  % i
            t1 = data(idxs,6);  % j
            t2 = data(idxs,7);  % i
            t3 = data(idxs,8);  % j
            t4 = data(idxs,9);  % i
            t5 = data(idxs,10); % j
            % calculate processed measurements
            T_rnd0 = t3-t0;
            T_rnd1 = t5-t2;
            T_rsp0 = t2-t1;
            T_rsp1 = t4-t3;
            %big Rij
            r = (obj.LIGHTSPEED)*(T_rnd0.*T_rnd1 - T_rsp0.*T_rsp1)./(T_rnd0 + T_rnd1 + T_rsp0 + T_rsp1);
            %r = (obj.LIGHTSPEED/2)*(T_rnd1 - T_rsp1);
            r1 =[];
            % rij if needed
            %obj.r_ij = (obj.LIGHTSPEED/2)*(T_rnd1 - T_rsp1);
            for i=1:length(r)
               if ~isnan(r(i))
                   r1=[r1;r(i)];
               end
            end
        end
        
        % get signal quality in a given time period
        function [fploss] = getFplossInWindow(obj, srcId, dstId, tstart, tstop)
            data = obj.dataparser.timinglogs;
            idxs = find(data(:,2) == srcId & data(:,3) == dstId & data(:,1) >= tstart & data(:,1) <= tstop);
            fploss = data(idxs,11);
        end
        
        % =============================================
        %              POSITION ACCESSORS
        % =============================================
        function xyz = getTruePosition(obj, nodeId, tarray)
            xyz = [];
            for i=1:length(tarray)
                t = tarray(i);
            idx = obj.getNodeIdx( nodeId );
            if obj.nodes{idx}.isMobile()
                rbid = obj.nodes{idx}.getRigidBodyId();
                p = obj.dataparser.getMocapPos( rbid, t );
            else
                p = obj.nodes{ idx }.getTruePosition()';
            end
            xyz = [xyz; p];
            end
        end
        
        function setRigidBodyId(obj, nodeId, rbId)
            nodeIdx = obj.getNodeIdx( nodeId );
            obj.nodes{nodeIdx}.setRigidBodyId( rbId );
        end
        
    end
    
end

