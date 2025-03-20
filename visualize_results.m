function visualize_results(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, warmup_frames)
% VISUALIZE_RESULTS 可视化真实轨迹和各种估计轨迹
%   VISUALIZE_RESULTS(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, warmup_frames) 
%   绘制真实位置和各种算法估计位置随时间的变化，并标记预热阶段

% 如果未提供warmup_frames参数，默认为4
if nargin < 7
    warmup_frames = 4;
end

% 获取轨迹长度
n_frames = size(true_positions, 1);
time_axis = 1:n_frames;

% 创建图形
figure('Position', [100, 100, 1200, 800]);

% 1. 距离随时间变化
subplot(3, 2, 1);
plot(time_axis, true_positions(:, 1), 'k-', 'LineWidth', 2);
hold on;
plot(time_axis, rd_positions(:, 1), 'g:', 'LineWidth', 1.5);
plot(time_axis, combined_positions(:, 1), 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_positions(:, 1), 'b--', 'LineWidth', 1.5);
plot(time_axis, estimated_positions(:, 1), 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('距离 (m)');
title('距离随时间变化');
legend('真实值', '距离多普勒', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 2. 距离误差随时间变化
subplot(3, 2, 2);
rd_range_error = rd_positions(:, 1) - true_positions(:, 1);
combined_range_error = combined_positions(:, 1) - true_positions(:, 1);
omp_range_error = omp_positions(:, 1) - true_positions(:, 1);
kf_range_error = estimated_positions(:, 1) - true_positions(:, 1);

plot(time_axis, rd_range_error, 'g:', 'LineWidth', 1.5);
hold on;
plot(time_axis, combined_range_error, 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_range_error, 'b--', 'LineWidth', 1.5);
plot(time_axis, kf_range_error, 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('距离误差 (m)');
title('距离估计误差');
legend('距离多普勒', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 3. 方位角随时间变化
subplot(3, 2, 3);
plot(time_axis, true_positions(:, 2), 'k-', 'LineWidth', 2);
hold on;
plot(time_axis, music_positions(:, 2), 'g:', 'LineWidth', 1.5);
plot(time_axis, combined_positions(:, 2), 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_positions(:, 2), 'b--', 'LineWidth', 1.5);
plot(time_axis, estimated_positions(:, 2), 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('方位角 (度)');
title('方位角随时间变化');
legend('真实值', 'MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 4. 方位角误差随时间变化
subplot(3, 2, 4);
music_azimuth_error = music_positions(:, 2) - true_positions(:, 2);
combined_azimuth_error = combined_positions(:, 2) - true_positions(:, 2);
omp_azimuth_error = omp_positions(:, 2) - true_positions(:, 2);
kf_azimuth_error = estimated_positions(:, 2) - true_positions(:, 2);

plot(time_axis, music_azimuth_error, 'g:', 'LineWidth', 1.5);
hold on;
plot(time_axis, combined_azimuth_error, 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_azimuth_error, 'b--', 'LineWidth', 1.5);
plot(time_axis, kf_azimuth_error, 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('方位角误差 (度)');
title('方位角估计误差');
legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 5. 俯仰角随时间变化
subplot(3, 2, 5);
plot(time_axis, true_positions(:, 3), 'k-', 'LineWidth', 2);
hold on;
plot(time_axis, music_positions(:, 3), 'g:', 'LineWidth', 1.5);
plot(time_axis, combined_positions(:, 3), 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_positions(:, 3), 'b--', 'LineWidth', 1.5);
plot(time_axis, estimated_positions(:, 3), 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('俯仰角 (度)');
title('俯仰角随时间变化');
legend('真实值', 'MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 6. 俯仰角误差随时间变化
subplot(3, 2, 6);
music_elevation_error = music_positions(:, 3) - true_positions(:, 3);
combined_elevation_error = combined_positions(:, 3) - true_positions(:, 3);
omp_elevation_error = omp_positions(:, 3) - true_positions(:, 3);
kf_elevation_error = estimated_positions(:, 3) - true_positions(:, 3);

plot(time_axis, music_elevation_error, 'g:', 'LineWidth', 1.5);
hold on;
plot(time_axis, combined_elevation_error, 'm-.', 'LineWidth', 1.5);
plot(time_axis, omp_elevation_error, 'b--', 'LineWidth', 1.5);
plot(time_axis, kf_elevation_error, 'r-', 'LineWidth', 1.5);

% 标记预热阶段
xline(warmup_frames + 0.5, '--', '预热结束', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
fill([0 warmup_frames+0.5 warmup_frames+0.5 0], [min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

grid on;
xlabel('帧');
ylabel('俯仰角误差 (度)');
title('俯仰角估计误差');
legend('MUSIC', '组合估计', 'OMP', '卡尔曼', 'Location', 'best');

% 调整整体布局
sgtitle('感知辅助太赫兹波束对准：参数估计性能', 'FontSize', 16);

% 创建3D轨迹图
figure('Position', [100, 100, 900, 700]);

% 将极坐标转换为笛卡尔坐标
true_cart = polar_to_cartesian(true_positions);
combined_cart = polar_to_cartesian(combined_positions);
omp_cart = polar_to_cartesian(omp_positions);
est_cart = polar_to_cartesian(estimated_positions);

% 绘制真实轨迹和估计轨迹
plot3(true_cart(:, 1), true_cart(:, 2), true_cart(:, 3), 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot3(combined_cart(:, 1), combined_cart(:, 2), combined_cart(:, 3), 'm-.s', 'LineWidth', 1.5, 'MarkerSize', 4);
plot3(omp_cart(:, 1), omp_cart(:, 2), omp_cart(:, 3), 'b--d', 'LineWidth', 1.5, 'MarkerSize', 4);
plot3(est_cart(:, 1), est_cart(:, 2), est_cart(:, 3), 'r-*', 'LineWidth', 1.5, 'MarkerSize', 6);

% 标记发射端位置
plot3(0, 0, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % 发射端位置

% 标记预热阶段和正常阶段的边界
if warmup_frames > 0 && warmup_frames < n_frames
    % 预热阶段结束点
    plot3(true_cart(warmup_frames, 1), true_cart(warmup_frames, 2), true_cart(warmup_frames, 3), ...
        'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'y');
    text(true_cart(warmup_frames, 1), true_cart(warmup_frames, 2), true_cart(warmup_frames, 3), ...
        '预热结束', 'FontSize', 12);
end

grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('接收端轨迹 (3D视图)', 'FontSize', 14);
legend('真实轨迹', '组合估计', 'OMP估计', '卡尔曼滤波', '发射端', 'Location', 'best');
view(30, 20); % 设置视角

% 添加坐标轴相等比例
axis equal;

end

function cart_coords = polar_to_cartesian(polar_coords)
% 将极坐标 [R, theta, phi] 转换为笛卡尔坐标 [x, y, z]
% R: 距离，theta: 方位角(度)，phi: 俯仰角(度)

n_points = size(polar_coords, 1);
cart_coords = zeros(n_points, 3);

for i = 1:n_points
    R = polar_coords(i, 1);
    theta = polar_coords(i, 2) * pi/180; % 转为弧度
    phi = polar_coords(i, 3) * pi/180;   % 转为弧度
    
    % 球坐标转笛卡尔坐标
    x = R * cos(phi) * sin(theta);
    y = R * cos(phi) * cos(theta);
    z = R * sin(phi);
    
    cart_coords(i, :) = [x, y, z];
end

end