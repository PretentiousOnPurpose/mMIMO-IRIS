
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Author(s): Yashwanth R - yashwanthr@iisc.ac.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear
close all;

rng(108);

[version, executable, isloaded] = pyversion;
if ~isloaded
    pyversion /usr/bin/python3
    py.print() %weird bug where py isn't loaded in an external script
end

% Params:
SIM_MOD                 = 1;

keepTXRX_Running = 1;

fprintf("5G Testbed@ IISc: IRIS MIMO Setup\n");
% fprintf("Transmission type: %s \n",chan_type);
fprintf("\n-------IRIS SDR 030 ----------\n");

%Iris params:
TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
chan_type               = "iris";
USE_HUB                 = 0;
TX_FRQ                  = 3.5e9;
RX_FRQ                  = TX_FRQ;
TX_GN                   = 42;
TX_GN_ue                = 42;
RX_GN                   = 42;
SMPL_RT                 = 5e6;
N_FRM                   = 10;
bs_ids                   = string.empty();
ue_ids                  = string.empty();
ue_scheds               = string.empty();

bs_ids = ["RF3E000064"]; %, "RF3E000037", "RF3E000012"];
numBSAnts = 2;

ue_ids= ["RF3E000044"];
numUEAnts = 2;

N_BS_NODE               = length(bs_ids);           % Number of nodes/antennas at the BS
N_UE                    = length(ue_ids);           % Number of UE nodes

numAntsPerBSNode = ceil(numBSAnts / N_BS_NODE);
numAntsPerUENode = ceil(numUEAnts / N_UE);

fprintf("Num of nodes per BS Nodes: %u\n", N_BS_NODE);
fprintf("Num of nodes per UE: %u\n", N_UE);
fprintf("Num of Antennas per BS Node: %u\n", numAntsPerBSNode);
fprintf("Num of Antennas per UE Node: %u\n", numAntsPerUENode);

while (1)
    %% Uplink Channel Estimation - UE
    
    % System Parameters
    iris_pre_zeropad = 0;
    iris_post_zeropad = 0;
    
    N = 256;
    CP = N/4;
    
    prmSeq = [zadoffChuSeq(5, (N+CP)-1); 0];
%     prmSeq = [prmSeq; prmSeq];
    
    PILOTS_SYMS_IND = [1:11];
%     DATA_SYMS_IND = [2;3;4;5;7;8;9;10];

    zeroSubcarriers = [127;128;129;];
    nonZeroSubcarriers = setdiff(1:N, zeroSubcarriers);
    numDataCarriers = length(nonZeroSubcarriers);
    
    numSyms = floor(4096 / (N + CP));
        
    ul_tx_data_buff_1 = zeros(N, numSyms - 1);
    ul_tx_data_buff_2 = zeros(N, numSyms - 1);
    
    ul_tx_data_buff_1(nonZeroSubcarriers, 1:2:11) = 1;
    ul_tx_data_buff_2(nonZeroSubcarriers, 2:2:11) = 1;
    
    % OFDM
    ul_tx_data_buff_1 = circshift(ul_tx_data_buff_1, N/2);
    ul_tx_data_buff_1 = ifft(ul_tx_data_buff_1, N);
    ul_tx_data_buff_1 = [ul_tx_data_buff_1(end-CP+1: end, :); ul_tx_data_buff_1];

    ul_tx_data_buff_2 = circshift(ul_tx_data_buff_2, N/2);
    ul_tx_data_buff_2 = ifft(ul_tx_data_buff_2, N);
    ul_tx_data_buff_2 = [ul_tx_data_buff_2(end-CP+1: end, :); ul_tx_data_buff_2];

    ul_tx_data(1, :) = [zeros(1, iris_pre_zeropad); prmSeq(:); ul_tx_data_buff_1(:); zeros(1, iris_post_zeropad)];
    ul_tx_data(2, :) = [zeros(1, iris_pre_zeropad); prmSeq(:); ul_tx_data_buff_2(:); zeros(1, iris_post_zeropad)];    
    
    %% Uplink Channel Estimation - gNB

    if (SIM_MOD == 1)
        delF = 0.0;
%         h = [[1,0];[0,1]];
        h = (1/sqrt(2)) * (randn(2, 2) + 1i * randn(2, 2));
        n = (1/sqrt(2*3.16*N)) * (randn(size(ul_tx_data)) + 1i * randn(size(ul_tx_data)));
        cfo_phase = exp(1i * 2 * pi * (delF) * (0:length(ul_tx_data)-1)' ./ N);
        ul_rx_data = h  * ul_tx_data + n;    
    end
    
     % Primary Sync.
    ul_rx_prm_corr = abs(xcorr(ul_rx_data(1, :), prmSeq));

    [~, max_idx] = max(ul_rx_prm_corr);
    max_idx = mod(max_idx, length(ul_tx_data)) + 1;
    
    % CFO Estimation and Equalization
%     rx_prmSeq = ul_rx_data(:, max_idx: max_idx + N + CP, 1);
    
%     corrSeq = (rx_prmSeq(161:320) .* conj(prmSeq(161:320))) .* conj((rx_prmSeq(1:160) .* conj(prmSeq(1:160))));
%     ffo_est = mean(angle(corrSeq)) * (256) / (2 * pi * 159);
%     ffo_est = 0;
% 
%     ul_rx_data = exp(-1i * 2 * pi * ffo_est * (0:length(ul_rx_data)-1)' ./ N) .* ul_rx_data;
    
    % Realign with Frame Boundary
    
    ul_rx_data = [ul_rx_data(:, 0 + max_idx + length(prmSeq): end), zeros(2, 0)];
    
    ul_rx_data_buff = zeros(2, N, numSyms - 1);
    
    % Framing
    
    for iter_syms = 1: (numSyms - 1)
        ul_rx_data_buff(1, :, iter_syms) = ul_rx_data(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP), 1);
        ul_rx_data_buff(2, :, iter_syms) = ul_rx_data(2, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP), 1);
    end
    
    % OFDM Demodulation
    
    ul_rx_data_buff_1 = fft(ul_rx_data_buff(1, :, :), N);
    ul_rx_data_buff_1 = squeeze(circshift(ul_rx_data_buff_1, N/2));

    ul_rx_data_buff_2 = fft(ul_rx_data_buff(2, :, :), N);
    ul_rx_data_buff_2 = squeeze(circshift(ul_rx_data_buff_2, N/2));

%     
%     % Phase Error Estimation and Correction
%     
%     delay1 = ul_rx_data_buff_1(nonZeroSubcarriers(2:end), :) .* conj(ul_rx_data_buff_1(nonZeroSubcarriers(1:end-1), :));
%     delay2 = ul_rx_data_buff_2(nonZeroSubcarriers(2:end), :) .* conj(ul_rx_data_buff_2(nonZeroSubcarriers(1:end-1), :));
%     
%     delay1 = -mean(angle(delay1(:))) * N / (2 * pi);
%     delay2 = -mean(angle(delay2(:))) * N / (2 * pi);
%     
%     ul_rx_data_buff_1(nonZeroSubcarriers, :) = exp(1i * 2 * pi * (nonZeroSubcarriers .* delay1)' ./ N) .* ul_rx_data_buff_1(nonZeroSubcarriers, :);
%     ul_rx_data_buff_2(nonZeroSubcarriers, :) = exp(1i * 2 * pi * (nonZeroSubcarriers .* delay2)' ./ N) .* ul_rx_data_buff_2(nonZeroSubcarriers, :);
    
    % Channel Estimation and Equalization
    
    ul_ch_est_11 = mean(mean(ul_rx_data_buff_1(nonZeroSubcarriers, 1:2:11), 2));
    ul_ch_est_12 = mean(mean(ul_rx_data_buff_1(nonZeroSubcarriers, 2:2:11), 2));
    ul_ch_est_21 = mean(mean(ul_rx_data_buff_2(nonZeroSubcarriers, 1:2:11), 2));
    ul_ch_est_22 = mean(mean(ul_rx_data_buff_2(nonZeroSubcarriers, 2:2:11), 2));
    
    H_hat_ul = [[ul_ch_est_11, ul_ch_est_12]; [ul_ch_est_21, ul_ch_est_22]];
     
    %% SVD based Spatial Multiplexing
    
    [U, S, V] = svd(H_hat_ul);
    [U1, S1, V1] = svd(h);
    
    %% Downlink Baseband Signal Generation - gNB
    
    % System Parameters
    iris_pre_zeropad = 0;
    iris_post_zeropad = 0;
    
    N = 256;
    CP = N/4;
    
    prmSeq = [zadoffChuSeq(5, (N+CP)-1); 0];
%     prmSeq = [prmSeq; prmSeq];
    
    codeRate = 1/2;
    modOrder = 4;
    
    PILOTS_SYMS_IND = [1;6;11];
    DATA_SYMS_IND = [2;3;4;5;7;8;9;10];

    zeroSubcarriers = [1;2;3;127;128;129;254;255;256];
    nonZeroSubcarriers = setdiff(1:N, zeroSubcarriers);
    numDataCarriers = length(nonZeroSubcarriers);
    
    numBits = (1/codeRate) * (1/log2(modOrder)) * length(DATA_SYMS_IND) * length(nonZeroSubcarriers);    
    numSyms = floor(4096 / (N + CP));

    % Data Generation
    
    dataBits(:, 1) = randi([0, 0], numBits, 1);
    dataBits(:, 2) = randi([1, 1], numBits, 1);
    
    % Channel coding
%   
    codedBits = zeros((1/codeRate) * numBits, 2);
    
    for iter_pc = 1: floor(numBits/100)
        codedBits_tmp = nrPolarEncode(dataBits((iter_pc - 1) * 100 + 1: iter_pc * 100, 1), (1/codeRate) * 100);
        codedBits((iter_pc - 1) * 200 + 1: iter_pc * 200, 1) = nrRateMatchPolar(codedBits_tmp, 100, (1/codeRate)*100);
        codedBits_tmp = nrPolarEncode(dataBits((iter_pc - 1) * 100 + 1: iter_pc * 100, 2), (1/codeRate) * 100);
        codedBits((iter_pc - 1) * 200 + 1: iter_pc * 200, 2) = nrRateMatchPolar(codedBits_tmp, 100, (1/codeRate)*100);
    end
    
    codedBits_tmp = nrPolarEncode(dataBits(iter_pc * 100 + 1: iter_pc * 100 + 76, 1), (1/codeRate) * 76);
    codedBits(iter_pc * 200 + 1: iter_pc * 200 + 152, 1) = nrRateMatchPolar(codedBits_tmp, 76, (1/codeRate)*76);
    codedBits_tmp = nrPolarEncode(dataBits(iter_pc * 100 + 1: iter_pc * 100 + 76, 2), (1/codeRate) * 76);
    codedBits(iter_pc * 200 + 1: iter_pc * 200 + 152, 2) = nrRateMatchPolar(codedBits_tmp, 76, (1/codeRate)*76);
        
%     codedBits = [dataBits; dataBits];
    % QAM
    
    modData = qammod(codedBits, modOrder, 'gray', 'InputType', 'bit', 'UnitAveragePower', 1);

    % SVD based SM Transmit Precoding
    modData = V * conj(modData');
    
    % Framing
        
    dl_tx_data_buff1 = zeros(N, numSyms - 1);
    dl_tx_data_buff2 = zeros(N, numSyms - 1);
    
    iter_data_syms = 1;
    iter_pilot_syms = 1;
    for iter_syms = 1: (numSyms - 1)
        if (any(DATA_SYMS_IND == iter_syms))        
            dl_tx_data_buff1(nonZeroSubcarriers, iter_syms) = modData(1, (iter_data_syms - 1) * numDataCarriers + 1: iter_data_syms * numDataCarriers);
            dl_tx_data_buff2(nonZeroSubcarriers, iter_syms) = modData(2, (iter_data_syms - 1) * numDataCarriers + 1: iter_data_syms * numDataCarriers);
            iter_data_syms = iter_data_syms + 1;
        else
%             dl_tx_data_buff1(nonZeroSubcarriers, iter_syms) = 1;
%             dl_tx_data_buff2(nonZeroSubcarriers, iter_syms) = 1;
%             iter_pilot_syms = iter_pilot_syms + 1;
        end            
    end
    
    dl_tx_data_buff1(nonZeroSubcarriers, PILOTS_SYMS_IND([1, 3])) = 1;
    dl_tx_data_buff2(nonZeroSubcarriers, PILOTS_SYMS_IND([2])) = 1;
    
    % OFDM
%     t1 = dl_tx_data_buff;
    dl_tx_data_buff1 = circshift(dl_tx_data_buff1, N/2);
    dl_tx_data_buff1 = ifft(dl_tx_data_buff1, N);
    dl_tx_data_buff1 = [dl_tx_data_buff1(end-CP+1: end, :); dl_tx_data_buff1];
    
    dl_tx_data_buff2 = circshift(dl_tx_data_buff2, N/2);
    dl_tx_data_buff2 = ifft(dl_tx_data_buff2, N);
    dl_tx_data_buff2 = [dl_tx_data_buff2(end-CP+1: end, :); dl_tx_data_buff2];


    dl_tx_data(1, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 320), reshape(dl_tx_data_buff1(:), 1, numel(dl_tx_data_buff1)), zeros(1, iris_post_zeropad)];
    dl_tx_data(2, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 320), reshape(dl_tx_data_buff2(:), 1, numel(dl_tx_data_buff1)), zeros(1, iris_post_zeropad)];
    
    %% Downlink Baseband Signal Processing at IRIS UE Nodes

    if (SIM_MOD == 1)
        delF = 0.0;
%         h = (1/sqrt(2)) * (randn(1, 1) + 1i * randn(1, 1));
        n = (1/sqrt(2*3.16*N)) * (randn(size(dl_tx_data)) + 1i * randn(size(dl_tx_data)));
        cfo_phase = exp(1i * 2 * pi * (delF) * (0:length(dl_tx_data)-1)' ./ N);
        dl_rx_data = h * dl_tx_data + n;    
    end
    
    % Primary Sync.
    dl_rx_prm_corr = abs(xcorr(dl_rx_data(1, :), prmSeq));

    [~, max_idx] = max(dl_rx_prm_corr);
    max_idx = mod(max_idx, length(dl_tx_data)) + 1;
    
    % CFO Estimation and Equalization
%     rx_prmSeq = dl_rx_data(max_idx: max_idx + N + CP, 1);
%     
%     corrSeq = (rx_prmSeq(161:320) .* conj(prmSeq(161:320))) .* conj((rx_prmSeq(1:160) .* conj(prmSeq(1:160))));
%     ffo_est = mean(angle(corrSeq)) * (256) / (2 * pi * 159);
%     
%     dl_rx_data = exp(-1i * 2 * pi * ffo_est * (0:length(dl_rx_data)-1)' ./ N) .* dl_rx_data;
    
    % Realign with Frame Boundary
%     
    dl_rx_data = dl_rx_data(:, max_idx + length(prmSeq): end);    
    dl_rx_data_buff1 = zeros(N, numSyms - 1);
    dl_rx_data_buff2 = zeros(N, numSyms - 1);
    
    % Framing
    
    for iter_syms = 1: (numSyms - 1)
        dl_rx_data_buff1(:, iter_syms) = dl_rx_data(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
        dl_rx_data_buff2(:, iter_syms) = dl_rx_data(2, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
    end
    
    % OFDM Demodulation
    
    dl_rx_data_buff1 = fft(dl_rx_data_buff1, N);
    dl_rx_data_buff1 = circshift(dl_rx_data_buff1, N/2);

    dl_rx_data_buff2 = fft(dl_rx_data_buff2, N);
    dl_rx_data_buff2 = circshift(dl_rx_data_buff2, N/2);
    
    % Channel Estimation and Equalization
    
    dl_ch_est_11 = mean(mean(dl_rx_data_buff1(nonZeroSubcarriers, PILOTS_SYMS_IND([1, 3])), 2));
    dl_ch_est_12 = mean(mean(dl_rx_data_buff1(nonZeroSubcarriers, PILOTS_SYMS_IND([2])), 2));
    dl_ch_est_21 = mean(mean(dl_rx_data_buff2(nonZeroSubcarriers, PILOTS_SYMS_IND([1, 3])), 2));
    dl_ch_est_22 = mean(mean(dl_rx_data_buff2(nonZeroSubcarriers, PILOTS_SYMS_IND([2])), 2));
    
    H_hat_dl = [[dl_ch_est_11, dl_ch_est_12]; [dl_ch_est_21, dl_ch_est_22]];
    [U2, S2, V2] = svd(H_hat_dl);
    % Equalized Data Extraction
    
    data_syms_buff1 = dl_rx_data_buff1(nonZeroSubcarriers, DATA_SYMS_IND);
    data_syms_buff2 = dl_rx_data_buff2(nonZeroSubcarriers, DATA_SYMS_IND);
    
    data_syms(1, :) = data_syms_buff1(:);
    data_syms(2, :) = data_syms_buff2(:);
    
    data_syms = U2' * data_syms;
    
    % QAM Demodulation
    
    demodData = qamdemod(conj(data_syms'), modOrder, 'gray', 'OutputType', 'approxllr', 'UnitAveragePower', 1);
%     decodedBits = demodData;
    % Channel Decoding
    
    decodedBits = zeros(numBits, 2);
    
    for iter_pc = 1: floor(numBits/100)
        decodedBits_tmp = nrRateRecoverPolar(demodData((iter_pc - 1) * 200 + 1: iter_pc * 200, 1), 100, 256);
        decodedBits((iter_pc - 1) * 100 + 1: iter_pc * 100, 1) = nrPolarDecode(decodedBits_tmp, 100, (1/codeRate)*100, 8);

        decodedBits_tmp = nrRateRecoverPolar(demodData((iter_pc - 1) * 200 + 1: iter_pc * 200, 2), 100, 256);
        decodedBits((iter_pc - 1) * 100 + 1: iter_pc * 100, 2) = nrPolarDecode(decodedBits_tmp, 100, (1/codeRate)*100, 8);
    end
    
    decodedBits_tmp = nrRateRecoverPolar(demodData(iter_pc * 200 + 1: iter_pc * 200 + 152, 1), 76, 256);
    decodedBits(iter_pc * 100 + 1: iter_pc * 100 + 76, 1) = nrPolarDecode(decodedBits_tmp, 76, (1/codeRate)*76, 8);
        
    decodedBits_tmp = nrRateRecoverPolar(demodData(iter_pc * 200 + 1: iter_pc * 200 + 152, 2), 76, 256);
    decodedBits(iter_pc * 100 + 1: iter_pc * 100 + 76, 2) = nrPolarDecode(decodedBits_tmp, 76, (1/codeRate)*76, 8);
    
    %% Evaluation 
%     fprintf("%u %u\n", sum(bitxor(dataBits(:, 1), decodedBits(:, 1))), sum(bitxor(dataBits(:, 2), decodedBits(:, 2))));
    if (sum(bitxor(dataBits, decodedBits)) == 0)
        fprintf("Successful Transmission\n");
    else
        fprintf("Failed Reception: %u\n", sum(bitxor(dataBits, decodedBits)));
    end
    
    if (keepTXRX_Running == 0)
        break;
    end
end


