function [range_est, velocity_est] = range_doppler_processing(rx_signal, params)
% RANGE_DOPPLER_PROCESSING 进行距离-多普勒处理和CFAR检测
%   [range_est, velocity_est] = RANGE_DOPPLER_PROCESSING(rx_signal, params)
%   对接收信号进行2D-FFT处理，生成距离-多普勒谱，并通过CFAR检测目标

% 获取信号维度
[n_rx_antennas, n_chirps, n_samples] = size(rx_signal);

% 初始化距离-多普勒谱
range_doppler_map = zeros(n_samples, n_chirps);

% 对子阵内所有天线的信号求平均 (非相干积累)
for ant_idx = 1:n_rx_antennas
    % 重新排列为[样本,chirp]格式
    ant_signal = squeeze(rx_signal(ant_idx, :, :))';
    
    % 应用窗函数 (减少旁瓣)
    range_window = hamming(n_samples);
    doppler_window = hamming(n_chirps);
    
    windowed_signal = ant_signal .* range_window;
    
    % 距离FFT (快时间FFT)
    range_fft = fft(windowed_signal, [], 1);
    
    % 对每个距离bin应用多普勒窗口
    for r_bin = 1:n_samples
        range_fft(r_bin, :) = range_fft(r_bin, :) .* doppler_window';
    end
    
    % 多普勒FFT (慢时间FFT)
    range_doppler = fft(range_fft, [], 2);
    
    % 非相干积累
    range_doppler_map = range_doppler_map + abs(range_doppler).^2;
end

% 平均化
range_doppler_map = range_doppler_map / n_rx_antennas;

% 调整频谱，使零多普勒位于中心
range_doppler_map = fftshift(range_doppler_map, 2);

% CFAR检测
[peaks, peak_indices] = cfar_detector(range_doppler_map, params);

if isempty(peaks)
    % 如果未检测到峰值，选择最强点
    [~, max_idx] = max(range_doppler_map(:));
    [r_idx, d_idx] = ind2sub(size(range_doppler_map), max_idx);
    peak_indices = [r_idx, d_idx];
end

% 取出最强峰值
[~, max_peak_idx] = max(peaks);
r_bin = peak_indices(max_peak_idx, 1);
d_bin = peak_indices(max_peak_idx, 2);

% 计算估计的距离
freq_res = params.fs / n_samples;
range_res = params.c / (2 * params.B);
if r_bin > n_samples/2
    r_bin = r_bin - n_samples; % 负频率校正
end
range_est = abs(r_bin) * range_res;

% 计算估计的速度
doppler_res = 1 / (n_chirps * params.T);
velocity_res = params.lambda / 2 * doppler_res;
d_bin_centered = d_bin - n_chirps/2 - 1; % 中心化
velocity_est = d_bin_centered * velocity_res;

end

function [peaks, peak_indices] = cfar_detector(range_doppler_map, params)
% 使用CA-CFAR (Cell-Averaging CFAR) 检测目标
% 参数解包
Pfa = params.Pfa;
guard_cells = params.guard_cells;
training_cells = params.training_cells;

% 获取数据尺寸
[n_range_bins, n_doppler_bins] = size(range_doppler_map);

% 初始化输出
peaks = [];
peak_indices = [];

% 计算CFAR窗口大小
range_guard = guard_cells(1);
doppler_guard = guard_cells(2);
range_train = training_cells(1);
doppler_train = training_cells(2);

% 设置从训练单元数量计算的CFAR缩放因子
% 根据恒虚警率的要求，我们计算缩放因子
num_training_cells = (2*range_train + 2*doppler_train + 4*range_train*doppler_train);
alpha = num_training_cells * (Pfa^(-1/num_training_cells) - 1);

% 对每个单元执行CFAR处理
for r_idx = 1+range_guard+range_train : n_range_bins-range_guard-range_train
    for d_idx = 1+doppler_guard+doppler_train : n_doppler_bins-doppler_guard-doppler_train
        % 当前单元值
        cell_under_test = range_doppler_map(r_idx, d_idx);
        
        % 提取训练单元区域 (排除保护单元)
        
        % 左侧训练区域
        left_train = range_doppler_map(r_idx-range_guard-range_train:r_idx-range_guard-1, ...
                                       d_idx-doppler_guard-doppler_train:d_idx+doppler_guard+doppler_train);
        
        % 右侧训练区域
        right_train = range_doppler_map(r_idx+range_guard+1:r_idx+range_guard+range_train, ...
                                        d_idx-doppler_guard-doppler_train:d_idx+doppler_guard+doppler_train);
        
        % 上方训练区域
        top_train = range_doppler_map(r_idx-range_guard:r_idx+range_guard, ...
                                      d_idx+doppler_guard+1:d_idx+doppler_guard+doppler_train);
        
        % 下方训练区域
        bottom_train = range_doppler_map(r_idx-range_guard:r_idx+range_guard, ...
                                         d_idx-doppler_guard-doppler_train:d_idx-doppler_guard-1);
        
        % 计算所有训练单元的平均值
        noise_level = mean([left_train(:); right_train(:); top_train(:); bottom_train(:)]);
        
        % 设置自适应阈值
        threshold = alpha * noise_level;
        
        % 检测是否超过阈值
        if cell_under_test > threshold
            % 局部最大值检测 (确保是峰值)
            local_region = range_doppler_map(max(1, r_idx-1):min(n_range_bins, r_idx+1), ...
                                             max(1, d_idx-1):min(n_doppler_bins, d_idx+1));
            if cell_under_test == max(local_region(:))
                peaks = [peaks; cell_under_test];
                peak_indices = [peak_indices; r_idx, d_idx];
            end
        end
    end
end

end 