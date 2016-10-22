%input_prefix = '../external_data/waiting_room/wear_20160929_171533';
input_prefix = '../external_data/waiting_room/wear_20160929_171859';

[acc, lacc, grav, gyro, mag, rot] = read_input(input_prefix);

clf

ax1 = subplot(4, 1, 1);
hold on
plot(acc(:, 1), acc(:, 2), 'r');
plot(acc(:, 1), acc(:, 3), 'g');
plot(acc(:, 1), acc(:, 4), 'b');
xlabel('acc')

ax2 = subplot(4, 1, 2);
hold on
plot(lacc(:, 1), lacc(:, 2), 'r');
plot(lacc(:, 1), lacc(:, 3), 'g');
plot(lacc(:, 1), lacc(:, 4), 'b');
xlabel('lacc')

ax3 = subplot(4, 1, 3);
hold on
plot(grav(:, 1), grav(:, 2), 'r');
plot(grav(:, 1), grav(:, 3), 'g');
plot(grav(:, 1), grav(:, 4), 'b');
xlabel('grav')

ax4 = subplot(4, 1, 4);
hold on
plot(gyro(:, 1), gyro(:, 2), 'r');
plot(gyro(:, 1), gyro(:, 3), 'g');
plot(gyro(:, 1), gyro(:, 4), 'b');
xlabel('gyro')

linkaxes([ax1, ax2, ax3, ax4], 'x')