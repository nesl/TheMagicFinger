%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160929_171533';  % Bo pointing
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160929_171859';  % Bo noise
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160929_180724';  % Bo pointing2
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160930_110758';  % Amr pointing1
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160930_112841';  % Bo pointing3
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160930_113349';  % Debbie
%input_prefix = '../external_data/pointing_segmentation/pilot2/wear_20160930_114004';  % Debbie
input_prefix = '../external_data/pointing_segmentation/pilot3/bo3/wear_20161006_164749';  % Bo pointing4

[acc, lacc, grav, gyro, mag, rot] = read_input(input_prefix);

clf

lacc_mag = sum(lacc(:, 2:4) .^ 2, 2) .^ 0.5;
gyro_mag = sum(gyro(:, 2:4) .^ 2, 2) .^ 0.5;


ax1 = subplot(2, 1, 1);
hold on
plot(lacc(:, 1), lacc_mag, 'k');
plot(lacc(:, 1), lacc(:, 2), 'r');
plot(lacc(:, 1), lacc(:, 3), 'g');
plot(lacc(:, 1), lacc(:, 4), 'b');
xlabel('lacc')
grid on

ax2 = subplot(2, 1, 2);
hold on
plot(gyro(:, 1), gyro_mag, 'k');
xlabel('gyro')
grid on

BIN_WIDTH = 0.05;
et = max(lacc(end, 1), gyro(end, 1));
num_bins = floor(et / BIN_WIDTH);

LACC_MAG_THRES = 15;
GYRO_MAG_THRES = 6;

idx = lacc_mag > LACC_MAG_THRES;
bin_idx = floor(lacc(idx, 1) / BIN_WIDTH) + 1;
lacc_mask = false(1, num_bins);
lacc_mask(bin_idx) = true;

idx = gyro_mag > GYRO_MAG_THRES;
bin_idx = floor(gyro(idx, 1) / BIN_WIDTH) + 1;
gyro_mask = false(1, num_bins);
gyro_mask(bin_idx) = true;

survive = lacc_mask & gyro_mask;
hold on
x = find(survive) * BIN_WIDTH;
y = ones(size(x)) * 20;
plot(x, y, 'mo');

SENSITIVE_LACC_MAG_THRES = 5;

sensitive_sensing = false(1, num_bins);
for i = 1:num_bins
    st = (i-1) * BIN_WIDTH;
    et = i     * BIN_WIDTH;
    idx = st < lacc(:, 1) & lacc(:, 1) < et;
    if max(lacc_mag(idx)) > SENSITIVE_LACC_MAG_THRES
        sensitive_sensing(i) = true;
    end
end

x = find(sensitive_sensing) * BIN_WIDTH;
y = ones(size(x)) * 18;
plot(x, y, 'gs');

linkaxes([ax1, ax2], 'x')