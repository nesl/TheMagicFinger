input_prefix = '../external_data/data/pre_pilot/wear_20160909_121108';


[acc, grav, gyro, mag, rot] = read_input(input_prefix);

clf

ax1 = subplot(5, 1, 1);
hold on
plot(acc(:, 1), acc(:, 2), 'r');
plot(acc(:, 1), acc(:, 3), 'g');
plot(acc(:, 1), acc(:, 4), 'b');
xlabel('acc')

ax2 = subplot(5, 1, 2);
hold on
plot(grav(:, 1), grav(:, 2), 'r');
plot(grav(:, 1), grav(:, 3), 'g');
plot(grav(:, 1), grav(:, 4), 'b');
xlabel('grav')

ax3 = subplot(5, 1, 3);
hold on
plot(gyro(:, 1), gyro(:, 2), 'r');
plot(gyro(:, 1), gyro(:, 3), 'g');
plot(gyro(:, 1), gyro(:, 4), 'b');
xlabel('gyro')

ax4 = subplot(5, 1, 4);
hold on
plot(mag(:, 1), mag(:, 2), 'r');
plot(mag(:, 1), mag(:, 3), 'g');
plot(mag(:, 1), mag(:, 4), 'b');
xlabel('mag')

ax5 = subplot(5, 1, 5);
hold on
plot(rot(:, 1), rot(:, 5), 'm');
%plot(rot(:, 1), rot(:, 6), 'y');
plot(rot(:, 1), rot(:, 2), 'r');
plot(rot(:, 1), rot(:, 3), 'g');
plot(rot(:, 1), rot(:, 4), 'b');
xlabel('rot')

linkaxes([ax1, ax2, ax3, ax4, ax5], 'x')