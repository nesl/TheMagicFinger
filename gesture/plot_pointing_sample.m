input_prefix = '../external_data/sample_command/wear_20161009_202012';  % Bo pointing4

offset = 1;
window_size = 2.5;
padding = 0.0;

xrange = [0 window_size] - padding - 0.05;

[acc, lacc, grav, gyro, mag, rot] = read_input(input_prefix);

lacc_mag = sum(lacc(:, 2:4) .^ 2, 2) .^ 0.5;
gyro_mag = sum(gyro(:, 2:4) .^ 2, 2) .^ 0.5;

cfigure(14, 4);
plot(lacc(:, 1) - offset, lacc(:, 2), 'r', 'LineWidth', 2);
hold on
plot(lacc(:, 1) - offset, lacc(:, 3), 'g', 'LineWidth', 2);
plot(lacc(:, 1) - offset, lacc(:, 4), 'b', 'LineWidth', 2);
grid on
xlim(xrange)
ylim([-40 40])
legend({'x', 'y', 'z'}, 'Orientation', 'Horizontal', 'Fontsize', 14)
saveplot('out_figure/lacc_xyz')

cfigure(14, 4);
plot(lacc(:, 1) - offset, lacc_mag, 'k', 'LineWidth', 2);
grid on
xlim(xrange)
ylim([0 50])
legend({'mag'}, 'Orientation', 'Horizontal', 'Fontsize', 14)
saveplot('out_figure/lacc_mag')

cfigure(14, 4);
plot(gyro(:, 1) - offset, gyro(:, 2), 'r', 'LineWidth', 2);
hold on
plot(gyro(:, 1) - offset, gyro(:, 3), 'g', 'LineWidth', 2);
plot(gyro(:, 1) - offset, gyro(:, 4), 'b', 'LineWidth', 2);
grid on
xlim(xrange)
ylim([-15 15])
legend({'x', 'y', 'z'}, 'Orientation', 'Horizontal', 'Fontsize', 14)
saveplot('out_figure/gyro_xyz')

cfigure(14, 4);
plot(gyro(:, 1) - offset, gyro_mag, 'k', 'LineWidth', 2);
grid on
xlim(xrange)
ylim([0 15])
legend({'mag'}, 'Orientation', 'Horizontal', 'Fontsize', 14)
saveplot('out_figure/gyro_mag')