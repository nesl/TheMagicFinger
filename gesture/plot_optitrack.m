file_path = '../external_data/pointing_segmentation/Bharath1/mocap_2016-09-28 17:37:57.csv';

a = csvread(file_path);
a(:, 1) = a(:, 1) - a(1, 1);
a(:, 1)
figure
hold on
plot(a(:, 1), a(:, 3), 'r-x')
plot(a(:, 1), a(:, 4), 'g-x')
plot(a(:, 1), a(:, 5), 'b-x')
