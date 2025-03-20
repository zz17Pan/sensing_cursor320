%% 感知辅助太赫兹波束对准系统
% 混合球面波-平面波(HSPM)模型 + 距离/角度估计 + OMP稀疏重建 + 卡尔曼滤波
% 主脚本：配置参数并调用各个模块

clear;
close all;
clc;

%% 系统参数设置
% 信号参数
params.c = 3e8;                 % 光速(m/s)
params.fc = 77e9;               % 载波频率(Hz)
params.lambda = params.c/params.fc; % 波长(m)
params.B = 1e9;                 % 带宽(Hz)
params.T = 50e-6;               % 扫频时间(s)
params.mu = params.B/params.T;  % 调频斜率(Hz/s)
params.fs = 20e6;               % 采样频率(Hz)
params.N_samples = 512;         % 每个chirp的采样点数
params.N_chirps = 64;           % chirp数量
params.T_frame = params.N_chirps * params.T; % 帧时间(s)
params.dt = params.T_frame;     % 时间步长(与帧时间相同)

% 天线阵列配置
params.N_tx_subarrays = 4;      % 发射子阵总数量
params.N_rx_subarrays = 4;      % 接收子阵总数量
params.N_antennas_per_subarray = 16; % 每个子阵的天线数(4x4)
params.d = params.lambda/2;     % 天线间距
params.d_sub = 16*params.lambda; % 子阵间距

% 感知阶段配置
params.sensing_tx_subarray = 1;  % 用于感知的发射子阵索引
params.sensing_rx_subarray = 1;  % 用于感知的接收子阵索引

% 目标参数(初始值)
params.initial_R = 40;          % 初始距离(m)
params.initial_theta = 20;      % 初始方位角(度)
params.initial_phi = 5;         % 初始俯仰角(度)

% CFAR参数
params.Pfa = 1e-6;              % 虚警率
params.guard_cells = [4, 4];    % 保护单元[距离, 多普勒]
params.training_cells = [12, 12]; % 训练单元[距离, 多普勒]

% MUSIC算法参数
params.angle_grid_step = 2;     % 角度搜索步长(度)
params.theta_range = [-60, 60]; % 方位角搜索范围(度)
params.phi_range = [-30, 30];   % 俯仰角搜索范围(度)

% OMP参数
params.max_iterations = 5;      % 最大迭代次数
params.residual_thr = 0.1;      % 残差阈值
params.R_grid_step = 0.5;       % 距离网格步长(m)

% 仿真参数
params.n_frames = 30;           % 仿真帧数
params.snr_db = 10;             % 信噪比(dB)

% 初始化发射端和接收端天线阵列位置
tx_array = initialize_tx_array(params);
rx_array = initialize_rx_array(params);

%% 仿真主循环
% 初始化卡尔曼滤波器
kf = initialize_kalman_filter(params);

% 保存结果的数组
estimated_positions = zeros(params.n_frames, 3); % 估计的[R,theta,phi]
true_positions = zeros(params.n_frames, 3);      % 真实的[R,theta,phi]
omp_positions = zeros(params.n_frames, 3);       % OMP估计的[R,theta,phi]
rd_positions = zeros(params.n_frames, 3);        % 距离-多普勒估计的[R,0,0] - 只估计距离
music_positions = zeros(params.n_frames, 3);     % MUSIC估计的[0,theta,phi] - 只估计角度
combined_positions = zeros(params.n_frames, 3);  % 组合估计的[R,theta,phi]

% 定义预热帧数
warmup_frames = 5;  % 增加预热帧数到5帧，确保有足够的稳定性

% 定义平滑过渡参数
smooth_frames = 3;  % 预热后的平滑过渡帧数
smooth_weight = 0.8; % 平滑权重：较大权重使过渡更平缓

% 开始仿真
fprintf('===========================================\n');
fprintf('开始感知辅助太赫兹波束对准仿真...\n');
fprintf('使用发射子阵 %d 和接收子阵 %d 进行感知\n', ...
        params.sensing_tx_subarray, params.sensing_rx_subarray);
fprintf('前 %d 帧使用真实值作为状态估计\n', warmup_frames);
fprintf('后续 %d 帧使用平滑过渡\n', smooth_frames);
fprintf('===========================================\n');

for frame_idx = 1:params.n_frames
    fprintf('\n处理帧 %d/%d\n', frame_idx, params.n_frames);
    fprintf('-------------------------------------------\n');
    
    % 更新接收端位置(移动目标)
    rx_array = update_rx_position(rx_array, params, frame_idx);
    
    % 获取当前真实的距离和角度
    [true_R, true_theta, true_phi] = calculate_true_params(tx_array, rx_array);
    true_positions(frame_idx, :) = [true_R, true_theta, true_phi];
    
    % 生成发射信号
    tx_signal = generate_fmcw_signal(params);
    
    % 使用混合球面平面波模型生成接收信号
    rx_signal = simulate_hspm_channel(tx_signal, tx_array, rx_array, params);
    
    % 检查是否是预热帧
    if frame_idx <= warmup_frames
        fprintf('预热阶段 %d/%d: 使用真实参数作为估计值\n', frame_idx, warmup_frames);
        
        % 直接使用真实值作为估计结果
        R_est = true_R;
        theta_est = true_theta;
        phi_est = true_phi;
        
        % 记录各种估计结果
        rd_positions(frame_idx, :) = [true_R, 0, 0];         % 距离-多普勒(只有距离)
        music_positions(frame_idx, :) = [0, true_theta, true_phi]; % MUSIC(只有角度)
        combined_positions(frame_idx, :) = [true_R, true_theta, true_phi];
        omp_positions(frame_idx, :) = [true_R, true_theta, true_phi];
        
        % 更新卡尔曼滤波器 (使用真实值作为观测)
        z = [true_R; true_theta*pi/180; true_phi*pi/180]; % 真实值作为观测
        kf = kalman_filter_update(kf, z, params);
        
        % 从卡尔曼滤波器获取估计结果
        x_est = kf.x(1);
        y_est = kf.x(4);
        z_est = kf.x(7);
        
        % 计算估计的球坐标
        estimated_R = sqrt(x_est^2 + y_est^2 + z_est^2);
        estimated_theta = atan2(y_est, x_est) * 180/pi;
        estimated_phi = atan2(z_est, sqrt(x_est^2 + y_est^2)) * 180/pi;
        
        % 保存估计结果
        estimated_positions(frame_idx, :) = [estimated_R, estimated_theta, estimated_phi];
    else
        % 正常处理流程：从预热帧之后开始
        fprintf('正常处理阶段 - 帧 %d\n', frame_idx);
        
        % 2D-FFT处理和CFAR检测 (获取距离和速度估计)
        [range_est, velocity_est] = range_doppler_processing(rx_signal, params);
        rd_positions(frame_idx, :) = [range_est, 0, 0]; % 只保存距离值
        
        % MUSIC角度估计 (仅估计角度)
        [theta_est, phi_est] = music_angle_estimation(rx_signal, params);
        music_positions(frame_idx, :) = [0, theta_est, phi_est]; % 只保存角度值
        
        % 结合距离-多普勒和MUSIC结果
        combined_positions(frame_idx, :) = [range_est, theta_est, phi_est];
        
        % 显示初步估计结果
        fprintf('距离-多普勒估计: R=%.2fm, v=%.2fm/s\n', range_est, velocity_est);
        fprintf('MUSIC角度估计: θ=%.2f°, φ=%.2f°\n', theta_est, phi_est);
        
        % 处理平滑过渡期
        if frame_idx <= warmup_frames + smooth_frames
            % 计算当前平滑帧在过渡期的位置 (1到smooth_frames)
            smooth_idx = frame_idx - warmup_frames;
            fprintf('平滑过渡期 %d/%d\n', smooth_idx, smooth_frames);
            
            % 获取上一帧的真实值或组合估计值
            if frame_idx == warmup_frames + 1
                % 第一个过渡帧，使用最后一个预热帧的真实值
                prev_R = true_positions(frame_idx-1, 1);
                prev_theta = true_positions(frame_idx-1, 2);
                prev_phi = true_positions(frame_idx-1, 3);
            else
                % 使用前一帧的组合估计值
                prev_R = combined_positions(frame_idx-1, 1);
                prev_theta = combined_positions(frame_idx-1, 2);
                prev_phi = combined_positions(frame_idx-1, 3);
            end
            
            % 计算动态平滑权重 - 从高到低逐步过渡
            weight_prev = smooth_weight * (1 - (smooth_idx-1)/smooth_frames);
            weight_curr = 1.0 - weight_prev;
            
            fprintf('  平滑权重: 前一帧=%.2f, 当前帧=%.2f\n', weight_prev, weight_curr);
            
            % 使用加权平均进行平滑
            R_smooth = weight_prev * prev_R + weight_curr * range_est;
            theta_smooth = weight_prev * prev_theta + weight_curr * theta_est;
            phi_smooth = weight_prev * prev_phi + weight_curr * phi_est;
            
            % 更新组合估计结果
            combined_positions(frame_idx, :) = [R_smooth, theta_smooth, phi_smooth];
            
            fprintf('  平滑前估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', range_est, theta_est, phi_est);
            fprintf('  平滑后估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', R_smooth, theta_smooth, phi_smooth);
            
            % 使用平滑后的值作为稀疏重建的先验
            range_est = R_smooth;
            theta_est = theta_smooth;
            phi_est = phi_smooth;
        end
        
        % OMP稀疏重建(使用距离-多普勒和MUSIC结果作为先验)
        prior_params.R_est = range_est;
        prior_params.theta_est = theta_est;
        prior_params.phi_est = phi_est;
        prior_params.sigma_R = 1.0;        % 距离估计的标准差，减小提高精度
        prior_params.sigma_theta = 3.0;    % 方位角估计的标准差，减小范围
        prior_params.sigma_phi = 3.0;      % 俯仰角估计的标准差，减小范围
        
        [R_omp, theta_omp, phi_omp] = omp_sparse_reconstruction(rx_signal, prior_params, params);
        omp_positions(frame_idx, :) = [R_omp, theta_omp, phi_omp];
        
        % 卡尔曼滤波更新
        z = [R_omp; theta_omp*pi/180; phi_omp*pi/180]; % 观测向量
        kf = kalman_filter_update(kf, z, params);
        
        % 从卡尔曼滤波器获取笛卡尔坐标状态，并转换回球坐标
        x_est = kf.x(1);
        y_est = kf.x(4);
        z_est = kf.x(7);
        
        % 计算估计的球坐标 (R, theta, phi)
        estimated_R = sqrt(x_est^2 + y_est^2 + z_est^2);
        estimated_theta = atan2(y_est, x_est) * 180/pi;
        estimated_phi = atan2(z_est, sqrt(x_est^2 + y_est^2)) * 180/pi;
        
        % 保存估计结果
        estimated_positions(frame_idx, :) = [estimated_R, estimated_theta, estimated_phi];
    end
    
    % 显示当前帧的估计结果
    fprintf('-------------------------------------------\n');
    fprintf('真实参数: R=%.2fm, θ=%.2f°, φ=%.2f°\n', true_R, true_theta, true_phi);
    fprintf('距离多普勒(R)+MUSIC(θ,φ): R=%.2fm, θ=%.2f°, φ=%.2f°\n', ...
            combined_positions(frame_idx,1), combined_positions(frame_idx,2), combined_positions(frame_idx,3));
    fprintf('OMP估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', omp_positions(frame_idx,1), omp_positions(frame_idx,2), omp_positions(frame_idx,3));
    fprintf('卡尔曼滤波估计: R=%.2fm, θ=%.2f°, φ=%.2f°\n', estimated_R, estimated_theta, estimated_phi);
    fprintf('估计误差: ΔR=%.2fm, Δθ=%.2f°, Δφ=%.2f°\n', ...
            abs(estimated_R-true_R), abs(estimated_theta-true_theta), abs(estimated_phi-true_phi));
    
    % 显示目标速度估计
    vx_est = kf.x(2);
    vy_est = kf.x(5);
    vz_est = kf.x(8);
    v_est = sqrt(vx_est^2 + vy_est^2 + vz_est^2);
    fprintf('估计目标速度: v=%.2fm/s (vx=%.2f, vy=%.2f, vz=%.2f)\n', ...
            v_est, vx_est, vy_est, vz_est);
    fprintf('真实目标速度: v=%.2fm/s (vx=%.2f, vy=%.2f, vz=%.2f)\n', ...
            norm(rx_array.velocity), rx_array.velocity(1), rx_array.velocity(2), rx_array.velocity(3));
end

%% 性能评估和可视化
% 计算各种方法的RMSE (从第warmup_frames+smooth_frames+1帧开始，排除预热和过渡部分)
eval_frames = (warmup_frames + smooth_frames + 1):params.n_frames;

% 距离-多普勒+MUSIC组合估计
rmse_R_combined = sqrt(mean((true_positions(eval_frames,1) - combined_positions(eval_frames,1)).^2));
rmse_theta_combined = sqrt(mean((true_positions(eval_frames,2) - combined_positions(eval_frames,2)).^2));
rmse_phi_combined = sqrt(mean((true_positions(eval_frames,3) - combined_positions(eval_frames,3)).^2));

% OMP估计
rmse_R_omp = sqrt(mean((true_positions(eval_frames,1) - omp_positions(eval_frames,1)).^2));
rmse_theta_omp = sqrt(mean((true_positions(eval_frames,2) - omp_positions(eval_frames,2)).^2));
rmse_phi_omp = sqrt(mean((true_positions(eval_frames,3) - omp_positions(eval_frames,3)).^2));

% 卡尔曼滤波
rmse_R_kf = sqrt(mean((true_positions(eval_frames,1) - estimated_positions(eval_frames,1)).^2));
rmse_theta_kf = sqrt(mean((true_positions(eval_frames,2) - estimated_positions(eval_frames,2)).^2));
rmse_phi_kf = sqrt(mean((true_positions(eval_frames,3) - estimated_positions(eval_frames,3)).^2));

fprintf('\n===========================================\n');
fprintf('性能评估 (RMSE，不包括预热和过渡阶段):\n');
fprintf('-------------------------------------------\n');
fprintf('距离多普勒+MUSIC: 距离=%.3fm, 方位角=%.3f°, 俯仰角=%.3f°\n', ...
        rmse_R_combined, rmse_theta_combined, rmse_phi_combined);
fprintf('OMP: 距离=%.3fm, 方位角=%.3f°, 俯仰角=%.3f°\n', ...
        rmse_R_omp, rmse_theta_omp, rmse_phi_omp);
fprintf('卡尔曼滤波: 距离=%.3fm, 方位角=%.3f°, 俯仰角=%.3f°\n', ...
        rmse_R_kf, rmse_theta_kf, rmse_phi_kf);
fprintf('===========================================\n');

% 可视化结果 (跟踪性能)
visualize_results(true_positions, rd_positions, music_positions, combined_positions, omp_positions, estimated_positions, warmup_frames+smooth_frames); 
