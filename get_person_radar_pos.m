function [radar_pos, ground_truth] = get_person_radar_pos(person_pos_key, radar_pos_key)
switch person_pos_key
    case "1_0"
        ground_truth = [-3.6, 4.8];
    case "1_1"
        ground_truth = [-4.8, 4.8];
    case "1_2"
        ground_truth = [-6.0, 4.8];
    case "2_0"
        ground_truth = [-3.6, 3.6];
    case "2_1"
        ground_truth = [-4.8, 3.6];
    case "2_2"
        ground_truth = [-6.0, 3.6];
    otherwise
        error("wrong person position")
end

switch radar_pos_key
    case "1_0"
        radar_pos = [-0.6, 2.4];
    case "1-1"
        radar_pos = [-0.6, 1.2];
    case "1-2"
        radar_pos = [-0.6, 0.0];
    case "2-0"
        radar_pos = [0.0, 2.4];
    case "2-1"
        radar_pos = [0.0, 1.2];
    case "2-2"
        radar_pos = [0.0, 0.0];
    otherwise
        error("wrong radar position")
end
end