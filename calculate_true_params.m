function [R, theta, phi] = calculate_true_params(tx_array, rx_array)
% CALCULATE_TRUE_PARAMS 计算发射端和接收端之间的真实距离和角度
%   [R, theta, phi] = CALCULATE_TRUE_PARAMS(tx_array, rx_array) 
%   计算发射端中心到接收端中心的距离、方位角和俯仰角

% 获取发射端和接收端中心位置
tx_center = tx_array.position;
rx_center = rx_array.position;

% 计算相对位置向量 (从发射端指向接收端)
rel_pos = rx_center - tx_center;

% 计算直线距离
R = norm(rel_pos);

% 计算从发射端到接收端的方位角和俯仰角
% 将相对位置转换为球坐标 (r, azimuth, elevation)
% MATLAB的cart2sph函数返回的是(azimuth, elevation, r)
[azimuth, elevation, ~] = cart2sph(rel_pos(1), rel_pos(2), rel_pos(3));

% 转换为度数并保留原始方向
theta = azimuth * 180/pi;
phi = elevation * 180/pi;

% 输出详细的调试信息
fprintf('计算真实参数:\n');
fprintf('  发射端位置: [%.2f, %.2f, %.2f]m\n', tx_center(1), tx_center(2), tx_center(3));
fprintf('  接收端位置: [%.2f, %.2f, %.2f]m\n', rx_center(1), rx_center(2), rx_center(3));
fprintf('  相对位置向量: [%.2f, %.2f, %.2f]m\n', rel_pos(1), rel_pos(2), rel_pos(3));
fprintf('  距离: %.2fm, 方位角: %.2f°, 俯仰角: %.2f°\n', R, theta, phi);

end 