%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Yashwanth R - yashwanthr@iisc.ac.in         %
% Title    : 4x1 MISO Transmit Beamforming (Closed Loop) %
% Date     : 29-04-2021                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

% rng(111);
%    
     
% [version, executable, isloaded] = pyversion;
% py.importlib.import_module('iris_py')

% pyversion /usr/bin/python3
% if ~isloaded
%     py.print() %weird bug where py isn't loaded in an external script
% end

% Params:
SIM_MOD = 1;

SNRdB = 100000000000;
SNR = 10 ^ (SNRdB / 10);

keepTXRX_Running = 1;
iter = 1;
numOK = 0;
numNOK = 0;
fprintf("5G Testbed@ IISc: IRIS MIMO Setup\n");
% fprintf("Transmission type: %s \n",chan_type);
fprintf("\n--------IRIS SDR 030 ----------\n");

%Iris params:
TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
chan_type               = "iris";
USE_HUB                 = 0;
TX_FRQ                  = 2.5e9;
RX_FRQ                  = TX_FRQ;
TX_GN                   = 42;
TX_GN_ue                = 42;
RX_GN                   = 20;
SMPL_RT                 = 5e6;
N_FRM                   = 1;
bs_ids                   = string.empty();
ue_ids                  = string.empty();
ue_scheds               = string.empty();

bs_ids = ["RF3E000064"]; %, "RF3E000037", "RF3E000012"];
numBSAnts = 4;

ue_ids= ["RF3E000044"];
numUEAnts = 1;

N_BS_NODE               = length(bs_ids);           % Number of nodes/antennas at the BS
N_UE                    = length(ue_ids);           % Number of UE nodes

numAntsPerBSNode = ceil(numBSAnts / N_BS_NODE);
numAntsPerUENode = ceil(numUEAnts / N_UE);


fprintf("Num of nodes per BS Nodes: %u\n", N_BS_NODE);
fprintf("Num of nodes per UE: %u\n", N_UE);
fprintf("Num of Antennas per BS Node: %u\n", numAntsPerBSNode);
fprintf("Num of Antennas per UE Node: %u\n", numAntsPerUENode);

while (1)
    %% Uplink Channel Estimation
    % System Parameters
    iris_pre_zeropad = 250;
    iris_post_zeropad = 250;
   
    N = 256;
    CP = N/4;
   
    prmSeq = [zadoffChuSeq(5, (N+CP)/2-1); 0];

    PILOTS_SYMS_IND = 1:9;
    numSyms = length(PILOTS_SYMS_IND) + 1;
%     DATA_SYMS_IND = [2,3,4,5,6,7,8];

    zeroSubcarriers = [];
    nonZeroSubcarriers = setdiff(1:N, zeroSubcarriers);
   
    % Framing
    ul_tx_data_buff = zeros(N, numSyms - 1);
   
    iter_pilot_syms = 1;
    for iter_syms = 1: (numSyms - 1)
        ul_tx_data_buff(nonZeroSubcarriers, iter_syms) = 1;
        iter_pilot_syms = iter_pilot_syms + 1;
    end
      
    % OFDM
    ul_tx_data_buff = circshift(ul_tx_data_buff, N/2);
    ul_tx_data_buff = ifft(ul_tx_data_buff, N);
    ul_tx_data_buff = [ul_tx_data_buff(end-CP+1: end, :); ul_tx_data_buff];
   

    ul_tx_data = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(ul_tx_data_buff(:), 1, numel(ul_tx_data_buff)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];

    n_samp = length(ul_tx_data);

    %% Downlink Baseband Signal Processing at IRIS UE Nodes

    if (SIM_MOD == 1)
        h = (1/sqrt(2)) * (randn(4, 1) + 1i * randn(4, 1));        
        h = h ./ h;
        ul_rx_data = h * ul_tx_data;
        noise = (1/sqrt(2 * SNR * N)) * (randn(size(ul_rx_data)) + 1i * randn(size(ul_rx_data)));    
        ul_rx_data = ul_rx_data + noise;
    else
        bs_sched = ["BGGGGGRG"];           % BS schedule
        ue_sched = ["GGGGGGPG"];           % UE schedule        
       
        % Iris nodes' parameters
        bs_sdr_params = struct(...
            'id', bs_ids, ...
            'n_sdrs', N_BS_NODE, ...        % number of nodes chained together
            'txfreq', TX_FRQ, ...  
            'rxfreq', RX_FRQ, ...
            'txgain', TX_GN, ...
            'rxgain', RX_GN, ...
            'sample_rate', SMPL_RT, ...
            'n_samp', n_samp, ...          % number of samples per frame time.
            'n_frame', N_FRM, ...
            'tdd_sched', bs_sched, ...     % number of zero-paddes samples
            'n_zpad_samp', iris_pre_zeropad ...
            );

        ue_sdr_params = bs_sdr_params;
        ue_sdr_params.id =  ue_ids;
        ue_sdr_params.n_sdrs = N_UE;
        ue_sdr_params.txgain = TX_GN_ue;

        ue_sdr_params.tdd_sched = ue_sched;
       
        n_samp = bs_sdr_params.n_samp;
       
        node_bs = iris_py(bs_sdr_params,[]);        % initialize BS
        node_ue = iris_py(ue_sdr_params,[]);    % initialize UE

        node_ue.sdr_configgainctrl();

        node_bs.sdrsync();                 % synchronize delays only for BS

        node_bs.sdrrxsetup();
        node_ue.sdrrxsetup();             % set up reading stream

%         tdd_sched_index = 1; % for uplink only one frame schedule is sufficient
        node_bs.set_tddconfig(1, bs_sdr_params.tdd_sched); % configure the BS: schedule etc.
        node_ue.set_tddconfig(0, ue_sdr_params.tdd_sched);

        node_bs.sdr_setupbeacon();   % Burn beacon to the BS(1) RAM

        for i=1:N_BS_NODE
%             dl_tx_data(1, :) = dl_tx_data(1, :); % ./ max(abs(dl_tx_data(1, :)));
%             dl_tx_data(2, :) = dl_tx_data(2, :); % ./ max(abs(dl_tx_data(2, :)));
            dl_tx_data = 1 .* dl_tx_data;
            node_ue.sdrtx(dl_tx_data); %, numAntsPerBSNode);       % Burn data to the UE RAM
        end

        node_bs.sdr_activate_rx();          % activate reading stream
        node_ue.sdr_setcorr()              % activate correlator        
%         node_bs.sdrtrigger();           % set trigger to start the frame
       
        % Iris Rx

        [dl_rx_data, data0_len] = node_bs.sdrrx(n_samp); % read data

        node_bs.sdr_close();                % close streams and exit gracefully.
        node_ue.sdr_close();
        fprintf('Length of the received vector from HW: \tUE:%d\n', data0_len);    
    end
   
    % Primary Sync.
    ul_rx_prm_corr = abs(xcorr(ul_rx_data(1, :), prmSeq));
%     figure;
%     plot(dl_rx_prm_corr);
    [ul_max_corr, max_idx] = max(ul_rx_prm_corr);
    max_idx = mod(max_idx, (length(ul_tx_data))) + 1;

    % Realign with Frame Boundary
    ul_rx_data = ul_rx_data(:, max_idx + length(prmSeq): end);
   
%     if (length(dl_rx_data1) < 960 || length(dl_rx_data2) < 960)
%         dl_rx_data = [[dl_rx_data1, zeros(1, 960-length(dl_rx_data1))]; [dl_rx_data2, zeros(1, 960-length(dl_rx_data2))]];
%     end   
   
    ul_rx_data_buff = zeros(4, N, numSyms - 1);
   
    % Framing   
    for iter_syms = 1: (numSyms - 1)
        ul_rx_data_buff(1, :, iter_syms) = ul_rx_data(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
        ul_rx_data_buff(2, :, iter_syms) = ul_rx_data(2, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
        ul_rx_data_buff(3, :, iter_syms) = ul_rx_data(3, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
        ul_rx_data_buff(4, :, iter_syms) = ul_rx_data(4, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
    end
   
    % OFDM Demodulation
   
    h_hat = zeros(4, 1);
    
    for iter_ant = 1:4
        ul_rx_data_buff1 = squeeze(fft(ul_rx_data_buff(iter_ant, : ,:), N));
        ul_rx_data_buff1 = circshift(ul_rx_data_buff1, N/2);

        % Channel Estimation and Equalization
        h_hat(iter_ant, 1) = mean(mean(ul_rx_data_buff1(:, PILOTS_SYMS_IND)));    
    end
    %% Downlink Baseband Signal Generation - gNB
   
    % System Parameters
    iris_pre_zeropad = 250;
    iris_post_zeropad = 250;
   
    N = 256;
    CP = N/4;
   
    prmSeq = [zadoffChuSeq(5, (N+CP)/2-1); 0];
   
    codeRate = 1/2;
    modOrder = 2;
   
    PILOTS_SYMS_IND = [1;9];
    DATA_SYMS_IND = [2,3,4,5,6,7,8];

    zeroSubcarriers = [];
    nonZeroSubcarriers = setdiff(1:N, zeroSubcarriers);
    numDataCarriers = length(nonZeroSubcarriers);
   
    numBits = floor((codeRate) * (log2(modOrder)) * length(DATA_SYMS_IND) * length(nonZeroSubcarriers));    
    numSyms = length(union(PILOTS_SYMS_IND, DATA_SYMS_IND)) + 1;

    % Data Generation
   
    dataBits = randi([0, 1], numBits, 1);
   
    % Channel coding
    codedBits = zeros((1/codeRate) * numBits, 1);
   
    bitGroup = 64;
    lastBits = numBits - (floor(numBits/bitGroup) * bitGroup);
    for iter_pc = 1: floor(numBits/bitGroup)
        codedBits_tmp = nrPolarEncode(dataBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 1), (1/codeRate) * bitGroup);
        codedBits((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 1) = nrRateMatchPolar(codedBits_tmp, bitGroup, (1/codeRate)*bitGroup);
    end
    
    if (lastBits ~= 0)
        codedBits_tmp = nrPolarEncode(dataBits(iter_pc * bitGroup + 1: iter_pc * bitGroup + lastBits, 1), (1/codeRate) * lastBits);
        codedBits(iter_pc * 2 * bitGroup + 1: iter_pc * 2 * bitGroup + 2 * lastBits, 1) = nrRateMatchPolar(codedBits_tmp, bitGroup, (1/codeRate)*lastBits);
    end
    % QAM    
    modData = qammod(codedBits, modOrder, 'gray', 'InputType', 'bit', 'UnitAveragePower', 1);
    modData = modData.';
    % Framing
    dl_tx_data_buff = zeros(N, numSyms - 1);
   
    iter_data_syms = 1;
    iter_pilot_syms = 1;
    for iter_syms = 1: (numSyms - 1)
        if (any(DATA_SYMS_IND == iter_syms))        
            dl_tx_data_buff(nonZeroSubcarriers, iter_syms) = modData(1, (iter_data_syms - 1) * numDataCarriers + 1: iter_data_syms * numDataCarriers);
            iter_data_syms = iter_data_syms + 1;
        else
            dl_tx_data_buff(nonZeroSubcarriers, iter_syms) = 1;
            iter_pilot_syms = iter_pilot_syms + 1;
        end            
    end
      
    % OFDM
    t1 = dl_tx_data_buff;
    dl_tx_data_buff = circshift(dl_tx_data_buff, N/2);
    dl_tx_data_buff = ifft(dl_tx_data_buff, N);
    dl_tx_data_buff = [dl_tx_data_buff(end-CP+1: end, :); dl_tx_data_buff];
   

    dl_tx_data(1, :) = h_hat(1, 1) .* [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff(:), 1, numel(dl_tx_data_buff)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
    dl_tx_data(2, :) = h_hat(2, 1) .* [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff(:), 1, numel(dl_tx_data_buff)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
    dl_tx_data(3, :) = h_hat(3, 1) .* [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff(:), 1, numel(dl_tx_data_buff)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
    dl_tx_data(4, :) = h_hat(4, 1) .* [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff(:), 1, numel(dl_tx_data_buff)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];

    dl_tx_data = (dl_tx_data) ./ 2;
    
    n_samp = length(dl_tx_data);

    %% Downlink Baseband Signal Processing at IRIS UE Nodes

    if (SIM_MOD == 1)
%         h = (1/sqrt(2)) * (randn(1, 4) + 1i * randn(1, 4));        
        dl_rx_data = h.' * dl_tx_data;
        noise = (1/sqrt(2 * SNR * N)) * (randn(size(dl_rx_data)) + 1i * randn(size(dl_rx_data)));    
        dl_rx_data = dl_rx_data + noise;
    else
        bs_sched = ["BGGGGGRG"];           % BS schedule
        ue_sched = ["GGGGGGPG"];           % UE schedule        
       
        % Iris nodes' parameters
        bs_sdr_params = struct(...
            'id', bs_ids, ...
            'n_sdrs', N_BS_NODE, ...        % number of nodes chained together
            'txfreq', TX_FRQ, ...  
            'rxfreq', RX_FRQ, ...
            'txgain', TX_GN, ...
            'rxgain', RX_GN, ...
            'sample_rate', SMPL_RT, ...
            'n_samp', n_samp, ...          % number of samples per frame time.
            'n_frame', N_FRM, ...
            'tdd_sched', bs_sched, ...     % number of zero-paddes samples
            'n_zpad_samp', iris_pre_zeropad ...
            );

        ue_sdr_params = bs_sdr_params;
        ue_sdr_params.id =  ue_ids;
        ue_sdr_params.n_sdrs = N_UE;
        ue_sdr_params.txgain = TX_GN_ue;

        ue_sdr_params.tdd_sched = ue_sched;
       
        n_samp = bs_sdr_params.n_samp;
       
        node_bs = iris_py(bs_sdr_params,[]);        % initialize BS
        node_ue = iris_py(ue_sdr_params,[]);    % initialize UE

        node_ue.sdr_configgainctrl();

        node_bs.sdrsync();                 % synchronize delays only for BS

        node_bs.sdrrxsetup();
        node_ue.sdrrxsetup();             % set up reading stream

%         tdd_sched_index = 1; % for uplink only one frame schedule is sufficient
        node_bs.set_tddconfig(1, bs_sdr_params.tdd_sched); % configure the BS: schedule etc.
        node_ue.set_tddconfig(0, ue_sdr_params.tdd_sched);

        node_bs.sdr_setupbeacon();   % Burn beacon to the BS(1) RAM

        for i=1:N_BS_NODE
%             dl_tx_data(1, :) = dl_tx_data(1, :); % ./ max(abs(dl_tx_data(1, :)));
%             dl_tx_data(2, :) = dl_tx_data(2, :); % ./ max(abs(dl_tx_data(2, :)));
            dl_tx_data = 1 .* dl_tx_data;
            node_ue.sdrtx(dl_tx_data); %, numAntsPerBSNode);       % Burn data to the UE RAM
        end

        node_bs.sdr_activate_rx();          % activate reading stream
        node_ue.sdr_setcorr()              % activate correlator        
%         node_bs.sdrtrigger();           % set trigger to start the frame
       
        % Iris Rx

        [dl_rx_data, data0_len] = node_bs.sdrrx(n_samp); % read data

        node_bs.sdr_close();                % close streams and exit gracefully.
        node_ue.sdr_close();
        fprintf('Length of the received vector from HW: \tUE:%d\n', data0_len);    
    end
   
    % Primary Sync.
    dl_rx_prm_corr = abs(xcorr(dl_rx_data, prmSeq));
%     figure;
%     plot(dl_rx_prm_corr);
    [max_corr, max_idx] = max(dl_rx_prm_corr);
    max_idx = mod(max_idx, (length(dl_tx_data))) + 1;

    rx_prm_seq = conj(dl_rx_data(1, max_idx: max_idx + length(prmSeq) - 1)');

    % Realign with Frame Boundary
    dl_rx_data = dl_rx_data(1, max_idx + length(prmSeq): end);
   
%     if (length(dl_rx_data1) < 960 || length(dl_rx_data2) < 960)
%         dl_rx_data = [[dl_rx_data1, zeros(1, 960-length(dl_rx_data1))]; [dl_rx_data2, zeros(1, 960-length(dl_rx_data2))]];
%     end   
   
    dl_rx_data_buff = zeros(N, numSyms - 1);
   
    % Framing   
    for iter_syms = 1: (numSyms - 1)
        dl_rx_data_buff(:, iter_syms) = dl_rx_data(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
    end
   
    % OFDM Demodulation
   
    dl_rx_data_buff = fft(dl_rx_data_buff, N);
    dl_rx_data_buff = circshift(dl_rx_data_buff, N/2);
   
    % Channel Estimation and Equalization
    h_hat = mean(mean(dl_rx_data_buff(:, PILOTS_SYMS_IND)));

    dl_rx_data_buff = dl_rx_data_buff .* (conj(h_hat) / (abs(h_hat) ^ 2));

    % Equalized Data Extraction
   
    data_syms_buff = dl_rx_data_buff(nonZeroSubcarriers, DATA_SYMS_IND);
   
    data_syms(:, 1) = data_syms_buff(:);
   
    % QAM Demodulation   
    demodData = qamdemod(data_syms, modOrder, 'gray', 'OutputType', 'bit', 'UnitAveragePower', 1);

    decodedBits = zeros(numBits, 1);
   
    for iter_pc = 1: floor(numBits/bitGroup)
        decodedBits_tmp = nrRateRecoverPolar(demodData((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 1), bitGroup, 128);
        decodedBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 1) = nrPolarDecode(decodedBits_tmp, bitGroup, (1/codeRate)*bitGroup, 8);
    end
    if (lastBits ~= 0)
        decodedBits_tmp = nrRateRecoverPolar(demodData(iter_pc * 2 * bitGroup + 1: iter_pc * 2 * bitGroup + 2 * lastBits, 1), lastBits, 128);
        decodedBits(iter_pc * bitGroup + 1: iter_pc * bitGroup + lastBits, 1) = nrPolarDecode(decodedBits_tmp, bitGroup, (1/codeRate)*lastBits, 8);
    end
    %% Evaluation
%     fprintf("%u %u\n", sum(bitxor(dataBits(:, 1), decodedBits(:, 1))), sum(bitxor(dataBits(:, 2), decodedBits(:, 2))));
    if (sum(bitxor(dataBits(:), decodedBits(:))) == 0)
        fprintf("Trail: %u | 4x1 MISO Beamforming - Successful Transmission\n", iter);
        numOK = numOK + 1;
    else
        fprintf("Trail: %u | 4x1 MISO Beamforming - Failed Transmission: %f bits in error\n", iter, sum(bitxor(dataBits(:), decodedBits(:)))/(2 * numBits));
        numNOK = numNOK + 1;
    end
   
    if (keepTXRX_Running == 0)
        break;
    end
   
    iter = iter + 1;
    
    if (iter == 15)
        break;
    end
end

fprintf("----------Simulation Results----------\n");
fprintf("Block Error Rate: %f", (numNOK)/(numOK + numNOK));
