classdef Task < handle
    % Task serves as a data structure to present the attributes of a task.
    
    properties (SetAccess = public, GetAccess = public)
        pathname;
        
        contributor;
        wearing_wrist;
        gnd_time_secs = [];
        gnd_events = [];
        
        % An array of Meta objects
        meta_recs = [];
    end
    
    methods
        % CONSTRUCTOR
        function obj = Task()
        end
        
        function [acc, grav, gyro, mag, rot] = get_sensor_data(obj)
            [acc, grav, gyro, mag, rot] = read_input(obj.pathname);
        end
        
        %{
        function sec = get_sensor_data_length(obj)
            exts = {'acc', 'gyro', 'grav', 'mag'};
            path_template = [obj.pathname '_%s.csv'];
            exp_col = 5;
            sensor_time_idx = 2;
            
            sec = inf;
            for i = 1:numel(exts)
                target_path = sprintf(path_template, exts{i});
                [ret_code_h, line_h] = system(['head -n ',num2str(1), ' ', target_path]);
                [ret_code_t, line_t] = system(['tail -n ',num2str(1), ' ', target_path]);
                if ret_code_h == 0 && ret_code_t == 0
                    mat_h = str2num(line_h);
                    mat_t = str2num(line_t);
                    if all(size(mat_h) == [1 exp_col]) && all(size(mat_t) == [1 exp_col])
                        cur_sec_diff = (mat_t(sensor_time_idx) - mat_h(sensor_time_idx)) * 1e-9;
                        sec = min(sec, cur_sec_diff);
                    end
                end
            end
            if isinf(sec)
                error(['Empty file or file does not exist: ' obj.pathname]);
            end
        end
        %}
        
        function res = required_attributes_all_set(obj)
            res = ~isempty(obj.pathname) & ~isempty(obj.contributor) ...
                & ~isempty(obj.wearing_wrist);
        end
    end
end