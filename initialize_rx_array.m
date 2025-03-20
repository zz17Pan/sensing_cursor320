function rx_array = initialize_rx_array(params)
% INITIALIZE_RX_ARRAY 初始化接收端天线阵列
%   rx_array = INITIALIZE_RX_ARRAY(params) 根据参数初始化接收端天线阵列结构
%   包括整体位置、子阵位置和子阵内天线位置

% 初始化接收端数组结构
rx_array = struct();

% 接收端初始位置 - 确保不在坐标轴上
initial_position = [25.0, 30.0, params.initial_R]; % 初始位置明确偏离坐标轴
rx_array.position = initial_position;

% 接收端运动参数 - 合理的速度矢量 (不要太极端)
rx_array.velocity = [12.0, -10.0, 8.0]; % 具有三个方向的速度分量

% 根据初始位置计算角度
theta = atan2(initial_position(2), initial_position(1)) * 180/pi; % 方位角
phi = atan2(initial_position(3), sqrt(initial_position(1)^2 + initial_position(2)^2)) * 180/pi; % 俯仰角

rx_array.theta = theta;
rx_array.phi = phi;

fprintf('初始化接收阵列:\n');
fprintf('  初始位置: [%.2f, %.2f, %.2f]m\n', ...
    initial_position(1), initial_position(2), initial_position(3));
fprintf('  方位角: %.2f°, 俯仰角: %.2f°\n', rx_array.theta, rx_array.phi);
fprintf('  速度矢量: [%.2f, %.2f, %.2f]m/s (总速度: %.2fm/s)\n', ...
    rx_array.velocity(1), rx_array.velocity(2), rx_array.velocity(3), ...
    norm(rx_array.velocity));

% 初始化子阵
rx_array.subarrays = cell(1, params.N_rx_subarrays);

% 为每个子阵分配位置和天线
for k = 1:params.N_rx_subarrays
    % 子阵中心位置 (沿x轴线性排列)
    subarray_center = initial_position + [(k - 2.5) * params.d_sub, 0, 0];
    
    % 初始化子阵结构
    rx_array.subarrays{k} = struct();
    rx_array.subarrays{k}.center = subarray_center;
    
    % 子阵内天线位置 (4x4 URA)
    antenna_positions = zeros(4, 4, 3); % [x,y,z]坐标
    
    % 计算每个天线位置
    for nx = 1:4
        for nz = 1:4
            antenna_positions(nx, nz, :) = subarray_center + ...
                [(nx - 2.5) * params.d, 0, (nz - 2.5) * params.d];
        end
    end
    
    rx_array.subarrays{k}.antenna_positions = antenna_positions;
end

end 
