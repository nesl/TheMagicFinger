%% Housekeeping
clc; close all; clear all;
add_paths

%% configuration
WIN_SIZE_SEC = 1;
WIN_SLIDE_SEC = 0.5;

TRAINING_RATIO = 0.7;


manager = get_task_manager();

rnd_seed = 5389;
rng(rnd_seed);

%% parsing input files
labels = [];
matrix = [];

for task_idx = 1:numel(manager.tasks)
    cur_task = manager.tasks{task_idx};
    
    [acc, grav, gyro, ~, ~] = cur_task.get_sensor_data();
    acc(:, 5) = (acc(:, 2) .^ 2 + acc(:, 3) .^ 2 + acc(:, 4) .^ 2) .^ 0.5;
    grav(:, 5) = (grav(:, 2) .^ 2 + grav(:, 3) .^ 2 + grav(:, 4) .^ 2) .^ 0.5;
    gyro(:, 5) = (gyro(:, 2) .^ 2 + gyro(:, 3) .^ 2 + gyro(:, 4) .^ 2) .^ 0.5;
    
    fprintf('parse file %s ... ', cur_task.pathname);
    if numrow(acc) == 0 || numrow(grav) == 0 || numrow(gyro) == 0
        fprintf('skipped\n');
        continue
    end
    
    stop_sec = min([acc(end, 1), grav(end, 1), gyro(end, 1)]);
    num_windows = floor((stop_sec - WIN_SIZE_SEC) / WIN_SLIDE_SEC);
    
    cur_labels = ones(num_windows, 1) * C.G_NONE;
    cur_matrix = zeros(num_windows, 240);
    
    for e_idx = 1:numel(cur_task.gnd_time_secs)
        w_s_idx = floor(cur_task.gnd_time_secs(e_idx) / WIN_SLIDE_SEC) + 1;
        w_e_idx = min(num_windows, w_s_idx + 1);
        cur_labels(w_s_idx:w_e_idx) = cur_task.gnd_events(e_idx);
    end
    
    for r_idx = 1:num_windows
        st = (r_idx - 1) * WIN_SLIDE_SEC;
        et = st + WIN_SIZE_SEC;
        
        acc_idx = st <= acc(:, 1) & acc(:, 1) < et;
        acc_window = acc(acc_idx, 2:end);
        grav_idx = st <= grav(:, 1) & grav(:, 1) < et;
        grav_window = grav(grav_idx, 2:end);
        gyro_idx = st <= gyro(:, 1) & gyro(:, 1) < et;
        gyro_window = gyro(gyro_idx, 2:end);
        
        windows = {acc_window, grav_window, gyro_window};
        cnt = 1;
        for w_idx = 1:numel(windows)
            this_window = windows{w_idx};
            for col_idx = 1:numcol(this_window)
                this_col = this_window(:, col_idx);
                cur_matrix(r_idx, cnt) = mean(this_col);          cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = min(this_col);           cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = max(this_col);           cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = var(this_col);           cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = std(this_col);           cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = prctile(this_col, 10);   cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = prctile(this_col, 25);   cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = prctile(this_col, 50);   cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = prctile(this_col, 75);   cnt = cnt + 1;
                cur_matrix(r_idx, cnt) = prctile(this_col, 90);   cnt = cnt + 1;
                tmp_fft = abs(fft(this_col)) / numel(this_col);
                cur_matrix(r_idx, cnt:cnt+9) = tmp_fft(1:10);     cnt = cnt + 10;
            end
        end
    end
    labels = [labels; cur_labels];
    matrix = [matrix; cur_matrix];
    
    fprintf('done\n');
end

label_dist = histc(labels, 1:5);
max_rows_each_moving_label = floor(min(label_dist(1:4)) + std(label_dist(1:4)));
max_rows_none_label = max_rows_each_moving_label * 5;

num_rows_all = [min(label_dist(1:4), max_rows_each_moving_label); ...
                min(label_dist(5), max_rows_none_label)];
fprintf('Number of rows in each class (and rows to be selected)\n');
for i = 1:5
    fprintf('    %12s: %4d -> %-4d\n', C.get_gesture_movement_name(i), label_dist(i), num_rows_all(i));
end


%% training and testing
training_labels = [];
training_matrix = [];
testing_labels = [];
testing_matrix = [];

for class_idx = 1:5
    idx = find(labels == class_idx);
    idx_of_idx = randperm(numel(idx));
    idx = idx(idx_of_idx);
    idx = idx(1:num_rows_all(class_idx));
    m = round(numel(idx) * TRAINING_RATIO);
    training_idx = idx(1:m);
    testing_idx = idx(m+1:end);
    training_labels = [training_labels; labels(training_idx)];
    training_matrix = [training_matrix; matrix(training_idx, :)];
    testing_labels = [testing_labels; labels(testing_idx)];
    testing_matrix = [testing_matrix; matrix(testing_idx, :)];
end

model = svmtrain(training_labels, training_matrix, '-c 1 -g 0.07');
[predict_label, accuracy, dec_values] = svmpredict(testing_labels, testing_matrix, model)
