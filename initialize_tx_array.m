function tx_array = initialize_tx_array(params)
% INITIALIZE_TX_ARRAY 初始化发射端天线阵列
%   tx_array = INITIALIZE_TX_ARRAY(params) 根据参数初始化发射端天线阵列结构
%   包括子阵位置和子阵内天线位置

% 初始化发射端数组结构
tx_array = struct();

% 发射端固定在原点
tx_array.position = [0, 0, 0]; 

% 初始化子阵
tx_array.subarrays = cell(1, params.N_tx_subarrays);

% 为每个子阵分配位置和天线
for k = 1:params.N_tx_subarrays
    % 子阵中心位置 (沿x轴线性排列)
    subarray_center = [(k - 2.5) * params.d_sub, 0, 0];
    
    % 初始化子阵结构
    tx_array.subarrays{k} = struct();
    tx_array.subarrays{k}.center = subarray_center;
    
    % 子阵内天线位置 (4x4 URA)
    antenna_positions = zeros(4, 4, 3); % [x,y,z]坐标
    
    % 计算每个天线位置
    for nx = 1:4
        for nz = 1:4
            antenna_positions(nx, nz, :) = subarray_center + ...
                [(nx - 2.5) * params.d, 0, (nz - 2.5) * params.d];
        end
    end
    
    tx_array.subarrays{k}.antenna_positions = antenna_positions;
end

end 