function rx_array = initialize_rx_array(params)
% INITIALIZE_RX_ARRAY 初始化接收端天线阵列
%   rx_array = INITIALIZE_RX_ARRAY(params) 根据参数初始化接收端天线阵列结构
%   包括整体位置、子阵位置和子阵内天线位置

% 初始化接收端数组结构
rx_array = struct();

% 接收端初始位置 - 避免在坐标轴上
initial_position = [15.0, 20.0, params.initial_R]; % 初始位置不在坐标轴上
rx_array.position = initial_position;

% 接收端运动参数 - 三维较大速度分量
rx_array.velocity = [15.0, -23.0, 30.0]; % 增大速度向量 [vx, vy, vz] (m/s)

% 当前角度 - 根据初始位置计算，而非固定值
[azimuth, elevation, ~] = cart2sph(initial_position(1), initial_position(2), initial_position(3));
rx_array.theta = azimuth * 180/pi; % 方位角(度)
rx_array.phi = elevation * 180/pi;   % 俯仰角(度)

fprintf('初始化接收阵列 - 位置: [%.2f, %.2f, %.2f]m, 方位角: %.2f°, 俯仰角: %.2f°\n', ...
    initial_position(1), initial_position(2), initial_position(3), ...
    rx_array.theta, rx_array.phi);
fprintf('速度矢量: [%.2f, %.2f, %.2f]m/s\n', ...
    rx_array.velocity(1), rx_array.velocity(2), rx_array.velocity(3));

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