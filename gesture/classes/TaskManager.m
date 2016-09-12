classdef TaskManager < handle
    % TaskManager keeps a set of tasks.
    %
    % TaskManager support builder pattern, the expect usage is:
    %
    %     manager = TaskManager();
    %     manager.add_new_task()
    %         .set_task_path(...)
    %         .set_contributor(...)
    %         .set_wearing_wrist(...)
    %         .set_gnd_time_secs(...)
    %         .set_gnd_events(...)
    %         .commit_task();
    
    properties (SetAccess = private)
        % variable tasks is a 1D-array of Task objects.
        tasks = Task().empty(0, 0);
        
        tmp_task;
        
        %{
        % for indexing
        %
        % query_indexing is a map serve as a hash. the key is a
        % concatenation of contributor and type. the value is an Nx2
        % matrix, each row represents (task_idx, set_idx).
        query_indexing;
        %}
    end
    
    methods
        % CONSTRUCTOR
        function obj = TaskManager()
        end
        
        %
        % Task builder funcions
        %
        function obj = add_new_task(obj)
            obj.tmp_task = Task();
        end
        
        function obj = set_task_path(obj, path)
            obj.tmp_task.pathname = path;
        end
        
        function obj = set_contributor(obj, contributor)
            obj.tmp_task.contributor = contributor;
        end
        
        function obj = set_wearing_wrist(obj, wrist)
            obj.tmp_task.wearing_wrist = C.get_wearing_wrist_label(wrist);
        end
        
        function obj = set_gnd_time_secs(obj, secs)
            obj.tmp_task.gnd_time_secs = secs;
        end
        
        function obj = set_gnd_events(obj, events)
            tmp = zeros(1, numel(events));
            for i = 1:numel(events)
                tmp(i) = C.get_gesture_movement_label(events{i});
            end
            obj.tmp_task.gnd_events = tmp;
        end
        
        %{
        function obj = add_meta_on_did_finish_set_computation_in_autocorr(obj, set_start_sec, set_stop_sec, action)
            obj.tmp_task.add_meta_on_did_finish_set_computation_in_autocorr(set_start_sec, set_stop_sec, action);
        end
        %}
        
        function obj = commit_task(obj)
            if numel(obj.tmp_task) == 0
                error('Method add_new_task() should be called first...')
            elseif ~obj.tmp_task.required_attributes_all_set()
                error('Some of the required attributes are missing ...')
            else
                obj.tasks{end+1} = obj.tmp_task;
                obj.tmp_task = [];
            end
        end
        
        
        %
        % Filter functions
        %
        %{
        function keep_types(obj, kept_types, varargin)
            if any(strcmp(varargin, 'replace-with-dont-care'))
                flag_replace_with_dont_care = true;
            end
            
            tasks_to_keep = false(size(obj.tasks));
            for task_idx = 1:numel(obj.tasks)
                % if there's no any workout session, skip this one
                if numel(obj.tasks{task_idx}.sets) == 0
                    continue
                end
                
                task = obj.tasks{task_idx};
                set_to_keep = false(size(task.sets));
                for i = 1:numel(task.sets)
                    if any(task.sets(i).type == kept_types)
                        set_to_keep(i) = true;
                        continue;
                    end
                end
                
                if flag_replace_with_dont_care && any(set_to_keep)
                    for i = 1:numel(task.sets)
                        if set_to_keep(i) == false
                            task.sets(i).type = TYPE.dont_care;
                            set_to_keep(i) = true;
                        end
                    end
                end
                
                task.sets = task.sets(set_to_keep);
                if numel(task.sets) > 0
                    tasks_to_keep(task_idx) = true;
                end
            end

            obj.tasks = obj.tasks(tasks_to_keep);
        end
        %}
        
        % Args:
        %   kept_contributors has to be a cell array
        function keep_contributors(obj, kept_contributors)
            f = @(t) any(strcmp(t.contributor, kept_contributors));
            tasks_to_keep = arrayfun(f, obj.tasks);
            obj.tasks = obj.tasks(tasks_to_keep);
        end
        
        
        %
        % Querying
        %
        
        % When accessing the following query-series function, indexing()
        % has to be called to initialize internally
        %{
        function indexing(obj)
            obj.query_indexing = containers.Map('KeyType','char','ValueType','any');
            for task_idx = 1:numel(obj.tasks)
                cur_task = obj.tasks{task_idx};
                for set_idx = 1:numel(cur_task.sets)
                    key = [cur_task.contributor '|' num2str(cur_task.sets(set_idx).type)];
                    if ~isKey(obj.query_indexing, key)
                        obj.query_indexing(key) = [task_idx, set_idx];
                    else
                        obj.query_indexing(key) = ...
                            [obj.query_indexing(key); task_idx, set_idx];
                    end
                end
            end
        end
        %}
        
        %
        % get contributor list
        %
        
        % Returns:
        %   A cell of contributors
        function contributors = get_all_contributors(obj)
            contributors = {};
            for task_idx = 1:numel(obj.tasks)
                if ~any(strcmp(obj.tasks{task_idx}.contributor, contributors))
                    contributors{end+1} = obj.tasks{task_idx}.contributor;
                end
            end
        end
        
        
        %
        % Statistics report functions
        %
        %{
        function report_by_workout_type(obj)
            total_num_sets = 0;
            total_num_reps = 0;
            total_workout_secs = 0;
            has_type = [];
            
            for j = 1:numel(obj.tasks)
                sets = obj.tasks{j}.sets;
                for k = 1:numel(sets)
                    has_type(sets(k).type) = true;
                end
            end
                
            for type = 1:numel(has_type)
                if ~has_type(type)
                    continue
                end
                if type == TYPE.dont_care
                    continue
                end
                
                num_sets = 0;
                num_reps = 0;
                workout_secs = 0;
                for j = 1:numel(obj.tasks)
                    sets = obj.tasks{j}.sets;
                    for k = 1:numel(sets)
                        if sets(k).type == type
                            num_sets = num_sets + 1;
                            num_reps = num_reps + sets(k).num_reps;
                            workout_secs = workout_secs + sets(k).stop_sec - sets(k).start_sec;
                        end
                    end
                end
                
                workout_name = TYPE.get_name(type);
                
                fprintf('%s: %d sets, %d reps, %.1f secs in workout\n', ...
                    workout_name, num_sets, num_reps, workout_secs);
                
                total_num_sets = total_num_sets + num_sets;
                total_num_reps = total_num_reps + num_reps;
                total_workout_secs = total_workout_secs + workout_secs;
            end
            
            fprintf('Totally %d sets.\n', total_num_sets);
            fprintf('Totally %d reps.\n', total_num_reps);
            fprintf('Totally %.1f mins in workout.\n', total_workout_secs / 60);
        end
        %}
        
        %{
        function report_total_length(obj)
            total_secs = 0;
            
            fprintf('Parsing file ');
            last_str = '';
            for i = 1:numel(obj.tasks)
                %[repmat('\b', 1, numel(last_str)) obj.tasks{i}.pathname]
                fprintf([repmat('\b', 1, numel(last_str)) obj.tasks{i}.pathname]);
                last_str = obj.tasks{i}.pathname;
                %obj.tasks{i}.pathname
                total_secs = total_secs + obj.tasks{i}.get_sensor_data_length();
                
            end
            fprintf('\n');
            fprintf('Totally %.1f mins in total.\n', total_secs / 60);
        end
        %}
        
        %{
        function report_contributing_status(obj)
            num_types = 0;
            map_type = containers.Map('KeyType','int32','ValueType','int32');
            type_list = {};
            num_contributors = 0;
            map_contributor = containers.Map('KeyType','char','ValueType','int32');
            contributor_list = {};
            
            % inquire the type set and the contributor set
            for task_idx = 1:numel(obj.tasks)
                cur_task = obj.tasks{task_idx};
                if ~isKey(map_contributor, cur_task.contributor)
                    num_contributors = num_contributors + 1;
                    map_contributor(cur_task.contributor) = num_contributors;
                    contributor_name = cur_task.contributor;
                    if length(contributor_name) > 6
                        contributor_name = [contributor_name(1:4) '..'];
                    end
                    contributor_list{num_contributors} = contributor_name;
                end
                
                for i = 1:numel(cur_task.sets)
                    if ~isKey(map_type, cur_task.sets(i).type)
                        num_types = num_types + 1;
                        map_type(cur_task.sets(i).type) = num_types;
                        type_list{num_types} = TYPE.get_name_short(cur_task.sets(i).type);
                    end
                end
            end
            
            % count number of sets and reps
            set_cnt = zeros(num_types, num_contributors);
            rep_cnt = zeros(num_types, num_contributors);
            
            for task_idx = 1:numel(obj.tasks)
                cur_task = obj.tasks{task_idx};
                contributor_idx = map_contributor(cur_task.contributor);
                
                for i = 1:numel(cur_task.sets)
                    type_idx = map_type(cur_task.sets(i).type);
                    set_cnt(type_idx, contributor_idx) ...
                        = set_cnt(type_idx, contributor_idx) + 1;
                    rep_cnt(type_idx, contributor_idx) ...
                        = rep_cnt(type_idx, contributor_idx) + cur_task.sets(i).num_reps;
                end
            end
            
            % print
            fprintf('[Set statistics]\n')
            fprintf('                      ');  % 22 spaces
            for i_contri = 1:numel(contributor_list)  % index of contributor
                fprintf('%7s', contributor_list{i_contri});
            end
            fprintf('\n');
            
            for i_type = 1:numel(type_list)
                fprintf('%20s  ', type_list{i_type});
                for j_contri = 1:numel(contributor_list)
                    fprintf('%7d', set_cnt(i_type, j_contri));
                end
                fprintf('%7d\n', sum(set_cnt(i_type, :)));  % sum
            end
            
            fprintf('                      ');  % 22 spaces
            for i_contri = 1:numel(contributor_list)  % index of contributor
                fprintf('%7d', sum(set_cnt(:, i_contri)));
            end
            fprintf('%7d\n', sum(sum(set_cnt)));
            
            
            fprintf('\n');
            
            fprintf('[Rep statistics]\n')
            fprintf('                      ');  % 22 spaces
            for i_contri = 1:numel(contributor_list)  % index of contributor
                fprintf('%7s', contributor_list{i_contri});
            end
            fprintf('\n');
            
            for i_type = 1:numel(type_list)
                fprintf('%20s  ', type_list{i_type});
                for j_contri = 1:numel(contributor_list)
                    fprintf('%7d', rep_cnt(i_type, j_contri));
                end
                fprintf('%7d\n', sum(rep_cnt(i_type, :)));  % sum
            end
            
            fprintf('                      ');  % 22 spaces
            for i_contri = 1:numel(contributor_list)  % index of contributor
                fprintf('%7d', sum(rep_cnt(:, i_contri)));
            end
            fprintf('%7d\n', sum(sum(rep_cnt)));
        end
        %}
        
    end
end