function [theta_est, phi_est] = music_angle_estimation(rx_signal, params)
% MUSIC_ANGLE_ESTIMATION 使用MUSIC算法估计目标的方位角和俯仰角
%   [theta_est, phi_est] = MUSIC_ANGLE_ESTIMATION(rx_signal, params)
%   对接收信号进行MUSIC算法处理，估计目标的方位角和俯仰角

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 对每个距离bin执行距离FFT
range_fft = zeros(n_rx_antennas, n_chirps, n_samples);
for ant_idx = 1:n_rx_antennas
    ant_signal = squeeze(rx_signal(ant_idx, :, :));
    range_fft(ant_idx, :, :) = fft(ant_signal, [], 2);
end

% 计算距离谱
range_spectrum = mean(abs(range_fft).^2, 2);
range_spectrum = squeeze(mean(range_spectrum, 1));

% 找到距离谱中的最强峰值对应的距离bin
[~, max_range_bin] = max(range_spectrum);

% 从该距离bin提取信号
signals_at_peak = squeeze(range_fft(:, :, max_range_bin));

% 估计协方差矩阵
R = signals_at_peak * signals_at_peak' / n_chirps;

% 对协方差矩阵进行特征分解
[V, D] = eig(R);
eigenvalues = diag(D);

% 按特征值降序排序
[eigenvalues, idx] = sort(eigenvalues, 'descend');
V = V(:, idx);

% 假设有1个信号源
n_sources = 1;

% 区分信号子空间和噪声子空间
noise_subspace = V(:, n_sources+1:end);

% 创建角度搜索网格
theta_range = params.theta_range;
phi_range = params.phi_range;
theta_grid = theta_range(1):params.angle_grid_step:theta_range(2);
phi_grid = phi_range(1):params.angle_grid_step:phi_range(2);

% 初始化MUSIC谱
n_theta = length(theta_grid);
n_phi = length(phi_grid);
music_spectrum = zeros(n_theta, n_phi);

% 计算MUSIC谱
for i = 1:n_theta
    theta = theta_grid(i) * pi/180;  % 转换为弧度
    
    for j = 1:n_phi
        phi = phi_grid(j) * pi/180;  % 转换为弧度
        
        % 计算阵列响应向量 (仅用于单个子阵)
        a = array_response_vector(theta, phi, params);
        
        % 计算MUSIC谱
        proj = noise_subspace' * a;
        music_spectrum(i, j) = 1 / (proj' * proj);
    end
end

% 找到MUSIC谱中的峰值
[~, max_idx] = max(music_spectrum(:));
[theta_idx, phi_idx] = ind2sub([n_theta, n_phi], max_idx);

% 转换为角度估计
theta_est = theta_grid(theta_idx);
phi_est = phi_grid(phi_idx);

end

function a = array_response_vector(theta, phi, params)
% 计算阵列响应向量 (仅用于单个子阵)
% theta: 方位角(弧度)
% phi: 俯仰角(弧度)

% 子阵内天线数量
Nx = 4;  % 子阵x方向天线数
Nz = 4;  % 子阵z方向天线数
n_antennas = Nx * Nz;

% 初始化阵列响应向量
a = zeros(n_antennas, 1);

% 计算方向余弦
u = cos(phi) * sin(theta);
v = sin(phi);
w = cos(phi) * cos(theta);

% 生成子阵内每个天线的响应
ant_idx = 1;
for nz = 1:Nz
    for nx = 1:Nx
        % 天线相对于子阵中心的位置
        x = (nx - 2.5) * params.d;
        z = (nz - 2.5) * params.d;
        
        % 计算相位
        phase = 2 * pi / params.lambda * (x * u + z * w);
        
        % 计算阵列响应
        a(ant_idx) = exp(1j * phase);
        
        ant_idx = ant_idx + 1;
    end
end

% 归一化
a = a / norm(a);

end 