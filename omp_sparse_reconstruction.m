function [R_est, theta_est, phi_est] = omp_sparse_reconstruction(rx_signal, prior_params, params)
% OMP_SPARSE_RECONSTRUCTION 使用先验信息的OMP算法进行参数估计
%   [R_est, theta_est, phi_est] = OMP_SPARSE_RECONSTRUCTION(rx_signal, prior_params, params)
%   利用距离-多普勒和MUSIC角度估计结果作为先验，通过OMP算法进行稀疏重建

% 获取先验信息
R_prior = prior_params.R_est;
theta_prior = prior_params.theta_est;
phi_prior = prior_params.phi_est;
sigma_R = prior_params.sigma_R;
sigma_theta = prior_params.sigma_theta;
sigma_phi = prior_params.sigma_phi;

% ===== 降低搜索范围和分辨率，减少计算量 =====
% 减小搜索范围 (以先验为中心，覆盖2个标准差而非3个)
R_range = [max(0.1, R_prior - 2*sigma_R), R_prior + 2*sigma_R];
theta_range = [theta_prior - 2*sigma_theta, theta_prior + 2*sigma_theta];
phi_range = [phi_prior - 2*sigma_phi, phi_prior + 2*sigma_phi];

% 确保角度范围在合理区间内
theta_range(1) = max(params.theta_range(1), theta_range(1));
theta_range(2) = min(params.theta_range(2), theta_range(2));
phi_range(1) = max(params.phi_range(1), phi_range(1));
phi_range(2) = min(params.phi_range(2), phi_range(2));

% 增大网格步长，减少字典尺寸
R_grid_step = params.R_grid_step * 2;    % 增大距离步长
angle_grid_step = params.angle_grid_step * 2; % 增大角度步长

% 定义参数网格
R_grid = R_range(1):R_grid_step:R_range(2);
theta_grid = theta_range(1):angle_grid_step:theta_range(2);
phi_grid = phi_range(1):angle_grid_step:phi_range(2);

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 降低信号采样维度 (可选)
downsample_factor = 2;
rx_signal_ds = rx_signal(:, 1:downsample_factor:end, 1:downsample_factor:end);
[~, n_chirps_ds, n_samples_ds] = size(rx_signal_ds);

% 重塑采样后的接收信号为矢量形式
y = reshape(rx_signal_ds, [], 1);

% 计算网格尺寸
n_R = length(R_grid);
n_theta = length(theta_grid);
n_phi = length(phi_grid);
n_atoms = n_R * n_theta * n_phi;

fprintf('字典参数: %d距离点 x %d方位角点 x %d俯仰角点 = %d个原子\n', ...
    n_R, n_theta, n_phi, n_atoms);

% ===== 内存高效的OMP实现 =====
% 不事先创建完整字典矩阵，而是按需计算

% 初始化
residual = y;
indices = [];
selected_params = [];
residual_norm = norm(residual);
initial_residual_norm = residual_norm;
iter = 0;

fprintf('开始OMP迭代，目标残差: %.2e\n', params.residual_thr * initial_residual_norm);

% OMP迭代
while (iter < params.max_iterations) && (residual_norm > params.residual_thr * initial_residual_norm)
    iter = iter + 1;
    fprintf('OMP迭代 %d/%d，当前残差: %.2e\n', iter, params.max_iterations, residual_norm);
    
    % 找到与残差最相关的原子
    max_correlation = -1;
    best_atom = [];
    best_atom_params = [];
    
    % 对每个可能的参数组合计算相关性
    for i_R = 1:n_R
        R = R_grid(i_R);
        fprintf('搜索距离 %.2f m (%d/%d)\n', R, i_R, n_R);
        
        for i_theta = 1:n_theta
            theta_deg = theta_grid(i_theta);
            theta = theta_deg * pi/180;
            
            for i_phi = 1:n_phi
                phi_deg = phi_grid(i_phi);
                phi = phi_deg * pi/180;
                
                % 计算权重 (高斯先验)
                w = exp(-(R - R_prior)^2/(2*sigma_R^2) - ...
                        (theta_deg - theta_prior)^2/(2*sigma_theta^2) - ...
                        (phi_deg - phi_prior)^2/(2*sigma_phi^2));
                
                % 生成字典原子
                atom = generate_atom(R, theta, phi, params, downsample_factor);
                
                % 应用先验权重
                atom = w * atom;
                
                % 计算与残差的相关性
                correlation = abs(atom' * residual);
                
                % 更新最大相关性
                if correlation > max_correlation
                    max_correlation = correlation;
                    best_atom = atom;
                    best_atom_params = [R, theta_deg, phi_deg];
                end
            end
        end
    end
    
    % 将最佳原子添加到选定集合
    indices = [indices; iter];
    selected_params = [selected_params; best_atom_params];
    
    if iter == 1
        % 第一次迭代
        Phi_selected = best_atom;
    else
        % 后续迭代
        Phi_selected = [Phi_selected, best_atom];
    end
    
    % 最小二乘估计
    s = Phi_selected \ y;
    
    % 更新残差
    residual = y - Phi_selected * s;
    residual_norm = norm(residual);
    
    fprintf('已选择参数: R=%.2f m, θ=%.2f°, φ=%.2f°\n', ...
        best_atom_params(1), best_atom_params(2), best_atom_params(3));
    fprintf('剩余残差: %.2e (%.2f%%)\n', residual_norm, 100*residual_norm/initial_residual_norm);
end

% 加权平均估计结果 (权重为稀疏系数的绝对值)
weights = abs(s);
weights = weights / sum(weights);

R_est = sum(weights .* selected_params(:, 1));
theta_est = sum(weights .* selected_params(:, 2));
phi_est = sum(weights .* selected_params(:, 3));

fprintf('OMP估计结果: R=%.2f m, θ=%.2f°, φ=%.2f°\n', R_est, theta_est, phi_est);
end

function atom = generate_atom(R, theta, phi, params, downsample_factor)
% 生成字典原子 (考虑降采样)
% R: 距离(m)
% theta: 方位角(弧度)
% phi: 俯仰角(弧度)
% downsample_factor: 降采样因子

% 使用单个接收子阵
n_rx_antennas = params.N_antennas_per_subarray;
n_chirps = params.N_chirps;
n_samples = params.N_samples;

% 降采样后的维度
n_chirps_ds = ceil(n_chirps / downsample_factor);
n_samples_ds = ceil(n_samples / downsample_factor);

% 计算原始分辨率下的时间点
t_chirp = (0:n_samples-1) / params.fs;

% 计算时延
tau = 2 * R / params.c;

% 方向矢量
direction = [cos(phi)*sin(theta), cos(phi)*cos(theta), sin(phi)];

% 估计多普勒频移
% 假设目标径向速度约为2m/s用于生成原子
v_r_assumed = 2.0;  
f_d = 2 * v_r_assumed * params.fc / params.c;

% 初始化原子
atom = zeros(n_rx_antennas, n_chirps_ds, n_samples_ds);

% 生成子阵内每个天线的响应
for ant_idx = 1:n_rx_antennas
    % 计算天线在子阵内的位置
    nx = mod(ant_idx-1, 4) + 1;
    nz = floor((ant_idx-1)/4) + 1;
    
    % 天线相对于子阵中心的位置
    x = (nx - 2.5) * params.d;
    z = (nz - 2.5) * params.d;
    
    % 空间相位
    phase = 2 * pi / params.lambda * (x*direction(1) + z*direction(3));
    
    % 计算接收信号 (考虑降采样)
    for chirp_idx = 1:downsample_factor:n_chirps
        ds_chirp_idx = ceil(chirp_idx / downsample_factor);
        if ds_chirp_idx > n_chirps_ds
            continue;
        end
        
        % 计算延迟后的时间
        t_delayed = t_chirp - tau;
        
        % 计算接收信号相位
        rx_phase = 2 * pi * (params.fc * t_delayed + 0.5 * params.mu * t_delayed.^2);
        
        % 添加多普勒频移和空间相位
        chirp_phase = rx_phase + 2*pi*f_d*t_chirp + phase;
        
        % 生成复指数信号 (降采样)
        for sample_idx = 1:downsample_factor:n_samples
            ds_sample_idx = ceil(sample_idx / downsample_factor);
            if ds_sample_idx > n_samples_ds
                continue;
            end
            atom(ant_idx, ds_chirp_idx, ds_sample_idx) = exp(1j * chirp_phase(sample_idx));
        end
    end
end

% 重塑为列向量
atom = reshape(atom, [], 1);

% 归一化
atom = atom / norm(atom);
end 