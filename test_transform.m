addpath(genpath('utils'));
theta = 45;
radar_pos = [-0.6, 2.4];
points = [0, 0; 0.6*sqrt(2), 0; 0, -0.6*sqrt(2)]; % 3x2
new_points = transform(points, radar_pos(1), radar_pos(2), 360-theta);
ground_truth = [-0.6, 2.4; 0, 3; 0, 1.8];
disp(new_points);
disp(ground_truth);
