%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Yashwanth R - yashwanthr@iisc.ac.in   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

% rng(111);
%    
     
[version, executable, isloaded] = pyversion;
py.importlib.import_module('iris_py')

% pyversion /usr/bin/python3
% if ~isloaded
%     py.print() %weird bug where py isn't loaded in an external script
% end

% Params:
SIM_MOD                 = 0;

keepTXRX_Running = 0;
iter = 1;
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
    %% Downlink Baseband Signal Generation - gNB
   
    % System Parameters
    iris_pre_zeropad = 300;
    iris_post_zeropad = 400;
   
    N = 256;
    CP = N/4;
   
    prmSeq = [zadoffChuSeq(5, (N+CP)/2-1); 0];
%     prmSeq = [prmSeq; prmSeq];
   
    codeRate = 1/2;
    modOrder = 4;
   
    PILOTS_SYMS_IND = [1;3];
    DATA_SYMS_IND = [2];

    zeroSubcarriers = [];
    nonZeroSubcarriers = setdiff(1:N, zeroSubcarriers);
    numDataCarriers = length(nonZeroSubcarriers);
   
    numBits = (1/codeRate) * (1/log2(modOrder)) * length(DATA_SYMS_IND) * length(nonZeroSubcarriers);    
    numSyms = length(union(PILOTS_SYMS_IND, DATA_SYMS_IND)) + 1;

    % Data Generation
   
    dataBits(:, 1) = randi([0, 1], numBits, 1);
    dataBits(:, 2) = randi([0, 1], numBits, 1);
   
    % Channel coding
    codedBits = zeros((1/codeRate) * numBits, 2);
   
    bitGroup = 64;
   
    for iter_pc = 1: floor(numBits/bitGroup)
        codedBits_tmp = nrPolarEncode(dataBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 1), (1/codeRate) * bitGroup);
        codedBits((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 1) = nrRateMatchPolar(codedBits_tmp, bitGroup, (1/codeRate)*bitGroup);
        codedBits_tmp = nrPolarEncode(dataBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 2), (1/codeRate) * bitGroup);
        codedBits((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 2) = nrRateMatchPolar(codedBits_tmp, bitGroup, (1/codeRate)*bitGroup);
    end


    % Convolutional Coding
%     codedBits = [dataBits(:, 1); dataBits(:, 1)];
%     codedBits(:, 2) = codedBits(:, 1);
%     trellis_a = struct('numInputSymbols',2,'numOutputSymbols',4, ...
% 'numStates',4,'nextStates',[0 2;0 2;1 3;1 3], ...
% 'outputs',[0 3;1 2;3 0;2 1]);
%     codedBits(:, 1) = convenc(dataBits(:, 1), trellis_a);
%     codedBits(:, 2) = convenc(dataBits(:, 2), trellis_a);

    % QAM    
    modData = qammod(codedBits, modOrder, 'gray', 'InputType', 'bit', 'UnitAveragePower', 1);

    % FFT based SM Transmit Precoding
    V = dftmtx(2);
    modData = V * conj(modData');
   
    % Framing
       
    dl_tx_data_buff1 = zeros(N, numSyms - 1);
%     dl_tx_data_buff2 = zeros(N, numSyms - 1);
   
    iter_data_syms = 1;
    iter_pilot_syms = 1;
    for iter_syms = 1: (numSyms - 1)
        if (any(DATA_SYMS_IND == iter_syms))        
            dl_tx_data_buff1(nonZeroSubcarriers, iter_syms) = modData(1, (iter_data_syms - 1) * numDataCarriers + 1: iter_data_syms * numDataCarriers);
            dl_tx_data_buff2(nonZeroSubcarriers, iter_syms) = modData(2, (iter_data_syms - 1) * numDataCarriers + 1: iter_data_syms * numDataCarriers);
            iter_data_syms = iter_data_syms + 1;
        else
            dl_tx_data_buff1(nonZeroSubcarriers, iter_syms) = 1;
            dl_tx_data_buff2(nonZeroSubcarriers, iter_syms) = 1;
            iter_pilot_syms = iter_pilot_syms + 1;
        end            
    end
   
%     dl_tx_data_buff1(nonZeroSubcarriers, PILOTS_SYMS_IND([1, 3])) = 1;
%     dl_tx_data_buff2(nonZeroSubcarriers, PILOTS_SYMS_IND([2, 4])) = 1;
   
    % OFDM
    t1 = dl_tx_data_buff1;
    dl_tx_data_buff1 = circshift(dl_tx_data_buff1, N/2);
    dl_tx_data_buff1 = ifft(dl_tx_data_buff1, N);
    dl_tx_data_buff1 = [dl_tx_data_buff1(end-CP+1: end, :); dl_tx_data_buff1];
   
    dl_tx_data_buff2 = circshift(dl_tx_data_buff2, N/2);
    dl_tx_data_buff2 = ifft(dl_tx_data_buff2, N);
    dl_tx_data_buff2 = [dl_tx_data_buff2(end-CP+1: end, :); dl_tx_data_buff2];


    dl_tx_data(1, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff1(:), 1, numel(dl_tx_data_buff1)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
    dl_tx_data(2, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(dl_tx_data_buff2(:), 1, numel(dl_tx_data_buff2)), 0.25 .* reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];

%     dl_tx_data(1, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
%     dl_tx_data(2, :) = [zeros(1, iris_pre_zeropad), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), reshape(prmSeq, 1, 160), zeros(1, iris_post_zeropad)];
%    

    n_samp = length(dl_tx_data);

    %% Downlink Baseband Signal Processing at IRIS UE Nodes

    if (SIM_MOD == 1)
        delF = 0.0;
        h = 1; %(1/sqrt(2)) * (randn(1, 1) + 1i * randn(1, 1));
        n = 0; %(1/sqrt(2*1*N)) * (randn(size(dl_tx_data)) + 1i * randn(size(dl_tx_data)));
        cfo_phase = exp(1i * 2 * pi * (delF) * (0:length(dl_tx_data)-1)' ./ N);
        dl_rx_data = h * dl_tx_data + n;    
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
    dl_rx_prm_corr1 = abs(xcorr(dl_rx_data(1, :), prmSeq));
    figure;
    subplot(1, 2, 1);
    plot(dl_rx_prm_corr1);
    [max_corr1, max_idx1] = max(dl_rx_prm_corr1);
    max_idx1 = mod(max_idx1, (length(dl_tx_data))) + 1;


    dl_rx_prm_corr2 = abs(xcorr(dl_rx_data(2, :), prmSeq));

    [max_corr2, max_idx2] = max(dl_rx_prm_corr2);
    max_idx2 = mod(max_idx2, length(dl_tx_data)) + 1;
    subplot(1, 2, 2);
    plot(dl_rx_prm_corr2);
 
    rx_prm_seq1 = conj(dl_rx_data(1, max_idx1: max_idx1 + length(prmSeq) - 1)');
    rx_prm_seq2 = conj(dl_rx_data(2, max_idx2: max_idx2 + length(prmSeq) - 1)');

    h1 = mean(rx_prm_seq1(1:end-1) ./ prmSeq(1:end-1));
    h2 = mean(rx_prm_seq2(1:end-1) ./ prmSeq(1:end-1));    
    % Realign with Frame Boundary
    dl_rx_data1 = dl_rx_data(1, max_idx1 + length(prmSeq): end); % .* (160 / max_corr1);
    dl_rx_data2 = dl_rx_data(2, max_idx2 + length(prmSeq): end); % .* (160 / max_corr2);
    dl_rx_data1 = dl_rx_data1 ./ h1;
    dl_rx_data2 = dl_rx_data2 ./ h2;
   
%     if (length(dl_rx_data1) < 960 || length(dl_rx_data2) < 960)
%         dl_rx_data = [[dl_rx_data1, zeros(1, 960-length(dl_rx_data1))]; [dl_rx_data2, zeros(1, 960-length(dl_rx_data2))]];
%     end
   
   
    dl_rx_data_buff1 = zeros(N, numSyms - 1);
    dl_rx_data_buff2 = zeros(N, numSyms - 1);
   
    % Framing
   
    for iter_syms = 1: (numSyms - 1)
        dl_rx_data_buff1(:, iter_syms) = dl_rx_data1(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
        dl_rx_data_buff2(:, iter_syms) = dl_rx_data2(1, (iter_syms - 1) * (N + CP) + CP + 1: iter_syms * (N + CP));
    end
   
    % OFDM Demodulation
   
    dl_rx_data_buff1 = fft(dl_rx_data_buff1, N);
    dl_rx_data_buff1 = circshift(dl_rx_data_buff1, N/2);

    dl_rx_data_buff2 = fft(dl_rx_data_buff2, N);
    dl_rx_data_buff2 = circshift(dl_rx_data_buff2, N/2);
   
    % Channel Estimation and Equalization

%     h1 = mean(mean(dl_rx_data_buff1(:, PILOTS_SYMS_IND)));
%     h2 = mean(mean(dl_rx_data_buff2(:, PILOTS_SYMS_IND)));
%    
%       dl_rx_data_buff1(:, DATA_SYMS_IND) = dl_rx_data_buff1(:, DATA_SYMS_IND) ./ h1;
%     dl_rx_data_buff1 = dl_rx_data_buff1 ./ h1;
%     dl_rx_data_buff2 = dl_rx_data_buff2 ./ h2;
%        
    % Equalized Data Extraction
   
    data_syms_buff1 = dl_rx_data_buff1(nonZeroSubcarriers, DATA_SYMS_IND);
    data_syms_buff2 = dl_rx_data_buff2(nonZeroSubcarriers, DATA_SYMS_IND);
   
    data_syms(1, :) = data_syms_buff1(:);
    data_syms(2, :) = data_syms_buff2(:);
   
%     data_syms = (data_syms(1, :) + data_syms(2, :)) ./ sqrt(2);
   
%     data_syms = (data_syms(1, :) + data_syms(2, :)) ./ (1.4142);
   
    % QAM Demodulation
   
    data_syms = V' * data_syms / 2;
   
    demodData = qamdemod(conj(data_syms'), modOrder, 'gray', 'OutputType', 'approxllr', 'UnitAveragePower', 1);

    decodedBits = zeros(numBits, 2);
   
    d1(:, 1) = dataBits(:, 1);
    d1(:, 2) = decodedBits(:, 1);
   
    d2(:, 1) = dataBits(:, 2);
    d2(:, 2) = decodedBits(:, 2);

   
    for iter_pc = 1: floor(numBits/bitGroup)
        decodedBits_tmp = nrRateRecoverPolar(demodData((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 1), bitGroup, 128);
        decodedBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 1) = nrPolarDecode(decodedBits_tmp, bitGroup, (1/codeRate)*bitGroup, 8);
        decodedBits_tmp = nrRateRecoverPolar(demodData((iter_pc - 1) * 2 * bitGroup + 1: iter_pc * 2 * bitGroup, 2), bitGroup, 128);
        decodedBits((iter_pc - 1) * bitGroup + 1: iter_pc * bitGroup, 2) = nrPolarDecode(decodedBits_tmp, bitGroup, (1/codeRate)*bitGroup, 8);
    end
    scatterplot(modData(:, 1));
    scatterplot(data_syms(:, 1));
    %% Evaluation
%     fprintf("%u %u\n", sum(bitxor(dataBits(:, 1), decodedBits(:, 1))), sum(bitxor(dataBits(:, 2), decodedBits(:, 2))));
    if (sum(bitxor(dataBits(:), decodedBits(:))) == 0)
        fprintf("Trail: %u | DFT-based Spatial Muliplexing - Successful Transmission\n", iter);
    else
        fprintf("Trail: %u | DFT-based Spatial Muliplexing - Failed Transmission: %f bits in error\n", iter, sum(bitxor(dataBits(:), decodedBits(:)))/(2 * numBits));
    end
   
    if (keepTXRX_Running == 0)
        break;
    end
   
    iter = iter + 1;
end
