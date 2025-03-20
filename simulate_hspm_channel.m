function rx_signal = simulate_hspm_channel(tx_signal, tx_array, rx_array, params)
% SIMULATE_HSPM_CHANNEL 使用简化的平面波模型模拟信道
%   rx_signal = SIMULATE_HSPM_CHANNEL(tx_signal, tx_array, rx_array, params)
%   仅使用指定的发射和接收子阵进行感知

% 获取信号形状
[n_chirps, n_samples] = size(tx_signal);

% 初始化接收信号(仅使用一个接收子阵的天线)
n_rx_antennas = params.N_antennas_per_subarray;
rx_signal = zeros(n_rx_antennas, n_chirps, n_samples);

% 计算每个chirp的时间点
t_chirp = (0:n_samples-1) / params.fs;

% 发射端和接收端之间的距离和角度
[R, theta_deg, phi_deg] = calculate_true_params(tx_array, rx_array);
theta = theta_deg * pi/180;
phi = phi_deg * pi/180;

% 计算径向速度(发射端到接收端方向的速度分量)
v_direction = [cos(phi)*sin(theta), cos(phi)*cos(theta), sin(phi)];
v_r = dot(rx_array.velocity, v_direction);

% 计算多普勒频移
f_d = 2 * v_r * params.fc / params.c;

% 获取用于感知的子阵
tx_subarray = tx_array.subarrays{params.sensing_tx_subarray};
rx_subarray = rx_array.subarrays{params.sensing_rx_subarray};

% 计算路径增益
path_gain = calculate_path_gain(R, params.fc, params.c);

% 生成接收子阵的方向矢量
a_rx = steeringVector2D(theta, phi, 4, 4, params.lambda, params.d);

% 处理每个chirp
for chirp_idx = 1:n_chirps
    % 将发射信号重塑为向量
    tx_vec = tx_signal(chirp_idx, :).';
    
    % 应用多普勒频移
    t = t_chirp;
    doppler_phase = exp(1j * 2 * pi * f_d * t);
    tx_vec_doppler = tx_vec .* doppler_phase.';
    
    % 通过信道传输信号
    for sample_idx = 1:n_samples
        % 计算时延对应的样本索引
        delay_samples = round(2 * R / params.c * params.fs);
        
        % 应用时延 (确保不超出边界)
        if sample_idx > delay_samples
            delayed_idx = sample_idx - delay_samples;
            % 使用方向矢量和路径增益生成接收信号
            rx_signal(:, chirp_idx, sample_idx) = path_gain * a_rx * tx_vec_doppler(delayed_idx);
        end
    end
end

% 添加噪声
snr_linear = 10^(params.snr_db/10);
signal_power = mean(abs(rx_signal(:)).^2);
noise_power = signal_power / snr_linear;
noise = sqrt(noise_power/2) * (randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
rx_signal = rx_signal + noise;

end

function path_gain = calculate_path_gain(distance, fc, c)
% 计算路径增益 (自由空间路径损耗模型)
lambda = c / fc;
path_gain = lambda / (4 * pi * distance) * exp(-1j * 2 * pi * distance / lambda);
end

function a = steeringVector2D(theta, phi, Nx, Nz, lambda, d)
% 生成二维子阵的方向矢量
a = zeros(Nx*Nz, 1);
idx = 1;
for nz = 0:(Nz-1)
    for nx = 0:(Nx-1)
        phase = 2*pi/lambda * (nx*d*sin(theta)*cos(phi) + nz*d*sin(phi));
        a(idx) = exp(1j*phase);
        idx = idx + 1;
    end
end
end 