
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Author(s): C. Nicolas Barati nicobarati@rice.edu 
%		Rahman Doost-Mohamamdy: doost@rice.edu
%
%
% Single-shot transmission from N_UE clients to N_BS_NODE base station
% radios (UE stands for User Equipment). We define two modes:
% OTA (Over-the-air) and SIM_MOD (simulation).
% In simulation mode we simply use a Rayleigh channel whereas the OTA mode
% relies on the Iris hardware for transmission and reception.
% In both cases the clients transmit an OFDM signal that resembles a
% typical 802.11 WLAN waveform. If the transmission is OTA, then the user
% specifies a schedule that tells all clients when to transmit their frame
% The base station initiates the schedule by sending a beacon signal that
% synchronizes clients. After that, all clients will transmit
% simultaneously. We implement a frame structure that allows the base
% station to capture clean (non-overlaping) training sequences for
% equalization and demultiplexing of the concurrent data streams.
%
%---------------------------------------------------------------------
% Original code copyright Mango Communications, Inc.
% Distributed under the WARP License http://warpproject.org/license
% Copyright (c) 2018-2019, Rice University
% RENEW OPEN SOURCE LICENSE: http://renew-wireless.org/license
% ---------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;

[version, executable, isloaded] = pyversion;
if ~isloaded
    pyversion /usr/bin/python
    py.print() %weird bug where py isn't loaded in an external script
end

% Params:
N_BS_NODE               = 1;
N_UE                    = 1;
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
SIM_MOD                 = 0;
DEBUG                   = 0;
PLOT                    = 1;
if SIM_MOD
    chan_type               = "rayleigh"; % Will use only Rayleigh for simulation
    sim_SNR_db              = 15;
    TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
    bs_ids                  = ones(1, N_BS_NODE);
    ue_ids                  = ones(1, N_UE);

else 
    %Iris params:
    TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
    chan_type               = "iris";
    USE_HUB                 = 0;
    TX_FRQ                  = 2.5e9;
    RX_FRQ                  = TX_FRQ;
    TX_GN                   = 42;
    TX_GN_ue                = 42;
    RX_GN                   = 42;
    SMPL_RT                 = 5e6;
    N_FRM                   = 10;
    bs_ids                   = string.empty();
    ue_ids                  = string.empty();
    ue_scheds               = string.empty();
    if USE_HUB
        % Using chains of different size requires some internal
        % calibration on the BS. This functionality will be added later.
        % For now, we use only the 4-node chains:

        bs_ids = ["RF3E000134", "RF3E000191", "RF3E000171", "RF3E000105",...
            "RF3E000053", "RF3E000177", "RF3E000192", "RF3E000117",...
            "RF3E000183", "RF3E000152", "RF3E000123", "RF3E000178", "RF3E000113", "RF3E000176", "RF3E000132", "RF3E000108", ...
            "RF3E000143", "RF3E000160", "RF3E000025", "RF3E000034",...
            "RF3E000189", "RF3E000024", "RF3E000139", "RF3E000032", "RF3E000154", "RF3E000182", "RF3E000038", "RF3E000137", ...
            "RF3E000103", "RF3E000180", "RF3E000181", "RF3E000188"];

        hub_id = "FH4A000001";

    else
        bs_ids = ["RF3E000008"]; %, "RF3E000017"]; %, "RF3E000024", "RF3E000139", "RF3E000032", "RF3E000154", "RF3E000182", "RF3E000038", "RF3E000137"];
    end

    ue_ids= ["RF3E000018"]; %, "RF3E000157"];

    N_BS_NODE               = length(bs_ids);           % Number of nodes/antennas at the BS
    N_UE                    = length(ue_ids);           % Number of UE nodes
end

fprintf("Channel type: %s \n",chan_type);


MIMO_ALG                = 'ZF';      % MIMO ALGORITHM: ZF or Conjugate 
fprintf("MIMO algorithm: %s \n",MIMO_ALG);

% Waveform params
N_OFDM_SYM              = 44;         % Number of OFDM symbols for burst, it needs to be less than 47
MOD_ORDER                        = 4;          % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)


% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
SC_IND_DATA_PILOT       = [2:27 39:64]';
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYM * length(SC_IND_DATA);       % Number of data symbols (one per data-bearing subcarrier per OFDM symbol) per UE
N_LTS_SYM               = 2;                                      % Number of 
N_SYM_SAMP              = N_SC + CP_LEN;                          % Number of samples that will go over the air
N_ZPAD_PRE              = 29;                                     % Zero-padding prefix for Iris
N_ZPAD_POST             = N_ZPAD_PRE -14;                         % Zero-padding postfix for Iris

% Rx processing params
FFT_OFFSET                    = 16;          % Number of CP samples to use in FFT (on average)
DO_APPLY_PHASE_ERR_CORRECTION = 0;           % Enable Residual CFO estimation/correction

%% Define the preamble
% LTS for fine CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 ...
    1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1].';
lts_t = ifft(lts_f, N_SC); %time domain

% Arrange time-orthogonal pilots 
preamble_common = [lts_t(33:64); repmat(lts_t,N_LTS_SYM,1)];
l_pre = length(preamble_common);
pre_z = zeros(size(preamble_common));
preamble = zeros(N_UE * l_pre, N_UE);
for jp = 1:N_UE
    preamble((jp-1)*l_pre + 1: (jp-1)*l_pre+l_pre,jp) = preamble_common;
end

%% Generate a payload of random integers
tx_data = randi(MOD_ORDER, N_DATA_SYMS, N_UE) - 1;

tx_syms = mod_sym(tx_data, MOD_ORDER);

% Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYM, N_UE);

% Define the pilot tone values as BPSK symbols
pilots = [1 1 -1 1].';
pilots = pilots .* 100;
% Repeat the pilots across all OFDM symbols
pilots_mat1 = repmat(pilots, 1, N_OFDM_SYM, N_UE);
pilots_mat2 = pilots_mat1;

pilots_mat1(:, 2:2:44) = 0;
pilots_mat2(:, 1:2:44) = 0;


%% IFFT

% Construct the IFFT input matrix
ifft_in_mat1 = zeros(N_SC, N_OFDM_SYM, N_UE);
ifft_in_mat2 = zeros(N_SC, N_OFDM_SYM, N_UE);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat1(SC_IND_DATA, :, :)   = tx_syms_mat;
ifft_in_mat1(SC_IND_PILOTS, :, :) = pilots_mat1;

ifft_in_mat2(SC_IND_DATA, :, :)   = tx_syms_mat;
ifft_in_mat2(SC_IND_PILOTS, :, :) = pilots_mat2;

%Perform the IFFT
tx_payload_mat1 = ifft(ifft_in_mat1, N_SC, 1);
tx_payload_mat2 = ifft(ifft_in_mat2, N_SC, 1);

% Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp1 = tx_payload_mat1((end-CP_LEN+1 : end), :, :);
    tx_payload_mat1 = [tx_cp1; tx_payload_mat1];

    tx_cp2 = tx_payload_mat2((end-CP_LEN+1 : end), :, :);
    tx_payload_mat2 = [tx_cp2; tx_payload_mat2];
end

% Reshape to a vector
tx_payload_vecs1 = reshape(tx_payload_mat1, ceil(numel(tx_payload_mat1)/N_UE), N_UE);
tx_payload_vecs2 = reshape(tx_payload_mat2, ceil(numel(tx_payload_mat2)/N_UE), N_UE);

% Construct the full time-domain OFDM waveform
tx_vecs1 = [zeros(N_ZPAD_PRE, N_UE); preamble; tx_payload_vecs1; zeros(N_ZPAD_POST, N_UE)];
tx_vecs2 = [zeros(N_ZPAD_PRE, N_UE); preamble; tx_payload_vecs2; zeros(N_ZPAD_POST, N_UE)];

% Leftover from zero padding:
tx_vecs_iris = [tx_vecs1, tx_vecs2];
% Scale the Tx vector to +/- 1
tx_vecs_iris = TX_SCALE .* tx_vecs_iris ./ max(abs(tx_vecs_iris));

%% SIMULATION:
if (SIM_MOD) 

    %number of samples in a frame
    n_samp = length(tx_vecs_iris);

    % Iris nodes' parameters
    bs_sdr_params = struct(...
        'id', bs_ids, ...
        'n_sdrs', N_BS_NODE, ...        % number of nodes chained together
        'txfreq', [], ...
        'rxfreq', [], ...
        'txgain', [], ...
        'rxgain', [], ...
        'sample_rate', [], ...
        'n_samp', n_samp, ...          % number of samples per frame time.
        'n_frame', [], ...
        'tdd_sched', [], ...     % number of zero-paddes samples
        'n_zpad_samp', N_ZPAD_PRE ...
        );
    
    ue_sdr_params = bs_sdr_params;
    ue_sdr_params.id =  ue_ids;
    ue_sdr_params.n_sdrs = N_UE;
    ue_sdr_params.txgain = [];

    rx_vec_iris = getRxVec(tx_vecs_iris, N_BS_NODE, N_UE, chan_type, sim_SNR_db, bs_sdr_params, ue_sdr_params, []);
    rx_vec_iris = rx_vec_iris.'; % just to agree with what the hardware spits out.
    
%% Init Iris nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the Iris experimenty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    % Create BS Hub and UE objects. Note: BS object is a collection of Iris
    % nodes.


    bs_sched = ["BGGGGGRG"];           % BS schedule
    ue_sched = ["GGGGGGPG"];           % UE schedule


    %number of samples in a frame
    n_samp = length(tx_vecs_iris);
   
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
        'n_zpad_samp', N_ZPAD_PRE ...
        );

    ue_sdr_params = bs_sdr_params;
    ue_sdr_params.id =  ue_ids;
    ue_sdr_params.n_sdrs = N_UE;
    ue_sdr_params.txgain = TX_GN_ue;
    
    ue_sdr_params.tdd_sched = ue_sched;
    
    if USE_HUB
        rx_vec_iris = getRxVec(tx_vecs_iris, N_BS_NODE, N_UE, chan_type, [], bs_sdr_params, ue_sdr_params, hub_id);
                               
    else
        rx_vec_iris = getRxVec(tx_vecs_iris, N_BS_NODE, N_UE, chan_type, [], bs_sdr_params, ue_sdr_params, []);
    end

end
% rx_vec_iris = sum(rx_vec_iris, 2) ./ sqrt(2);
l_rx_dec=length(rx_vec_iris);
% rx_vec_iris = rx_vec_iris(2, :);

N_BS_NODE = 2;
N_UE = 2;
%% Correlate for LTS
% Complex cross correlation of Rx waveform with time-domain LTS

a = 1;
unos = ones(size(preamble_common'));
lts_corr = zeros(N_BS_NODE, length(rx_vec_iris));
data_len = (N_OFDM_SYM)*(N_SC +CP_LEN);
rx_lts_mat = double.empty();
payload_ind = int32.empty();
payload_rx = zeros(N_BS_NODE, data_len);
lts_peaks = zeros(N_BS_NODE, N_UE);

for ibs =1:N_BS_NODE
        % Correlation through filtering
        v0 = filter(fliplr(preamble_common'),a,rx_vec_iris(ibs,:));
        v1 = filter(unos,a,abs(rx_vec_iris(ibs,:)).^2);
        lts_corr(ibs,:) = (abs(v0).^2)./v1; % normalized correlation
        
        % Sort the correlation values
        sort_corr = sort(lts_corr(ibs,:), 'descend');
        % Take the N_UE largest values
        rho_max = sort_corr(1:N_UE);
        % Get the indices of N_UE largest corr. values
        lts_peaks(ibs,:) = find(lts_corr(ibs,:) >= min(rho_max));
        
        % position of the last peak
        max_idx = max(lts_peaks(ibs,:));
        
        % In case of bad correlatons:
        if (max_idx + data_len) > length(rx_vec_iris) || (max_idx < 0) || (max_idx - length(preamble) < 0)
            fprintf('Bad correlation at antenna %d max_idx = %d \n', ibs, max_idx);
            % Real value doesn't matter since we have corrrupt data:
            max_idx = length(rx_vec_iris)-data_len -1;
            
        end
            
        payload_ind(ibs) = max_idx +1;
        pream_ind_ibs = payload_ind(ibs) - length(preamble);
        pl_idx = payload_ind(ibs) : payload_ind(ibs) + data_len;
        rx_lts_mat(ibs,:) = rx_vec_iris(ibs, pream_ind_ibs: pream_ind_ibs + length(preamble) -1 );
        payload_rx(ibs,1:length(pl_idx) -1) = rx_vec_iris(ibs, payload_ind(ibs) : payload_ind(ibs) + length(pl_idx) -2);

end

if DEBUG
    figure,
    for sp = 1:N_BS_NODE
        subplot(N_BS_NODE,1,sp);
        plot(lts_corr(sp,:));
        grid on;
        xlabel('Samples');
        y_label = sprintf('Anetnna %d',sp);
        ylabel(y_label);
    end
    sgtitle('LTS correlations accross antennas')
end



%% Rx processing

% Reshape the payload and take subcarriers without the CP
payload_rx = reshape(payload_rx,N_BS_NODE, (N_SC + CP_LEN),[]);
payload_noCP = payload_rx(:,CP_LEN-FFT_OFFSET+(1:N_SC),:);
% Take FFT
Y_data = fft(payload_noCP, N_SC, 2);
x = reshape(Y_data(1, :, :), 64, 44);

% Channel Estimation
pilots_11 = reshape(Y_data(1, SC_IND_PILOTS, (1:2:44)), 4, 22);
pilots_12 = reshape(Y_data(1, SC_IND_PILOTS, (2:2:44)), 4, 22);
pilots_21 = reshape(Y_data(2, SC_IND_PILOTS, (1:2:44)), 4, 22);
pilots_22 = reshape(Y_data(2, SC_IND_PILOTS, (2:2:44)), 4, 22);

pilots_11 = pilots_11 ./ pilots;
pilots_12 = pilots_12 ./ pilots;
pilots_21 = pilots_21 ./ pilots;
pilots_22 = pilots_22 ./ pilots;

% Estimating Timing Offset and correcting phase errors
phase_11 = mean(mean(mod(angle(pilots_11), 2 * pi), 2) ./ ((2 * pi / 64) * (SC_IND_PILOTS - 1)'), 1);
phase_12 = mean(mean(mod(angle(pilots_12), 2 * pi), 2) ./ ((2 * pi / 64) * (SC_IND_PILOTS - 1)'), 1);
phase_21 = mean(mean(mod(angle(pilots_21), 2 * pi), 2) ./ ((2 * pi / 64) * (SC_IND_PILOTS - 1)'), 1);
phase_22 = mean(mean(mod(angle(pilots_22), 2 * pi), 2) ./ ((2 * pi / 64) * (SC_IND_PILOTS - 1)'), 1);

D = mean([phase_11, phase_12, phase_21, phase_22]);

pilots_11 = mean(pilots_11, 2);
pilots_12 = mean(pilots_12, 2);
pilots_21 = mean(pilots_21, 2);
pilots_22 = mean(pilots_22, 2);

ch_11 = mean(pilots_11(:) .* exp(-1i * 2 * pi * (SC_IND_PILOTS - 1)' * D / 64));
ch_12 = mean(pilots_12(:) .* exp(-1i * 2 * pi * (SC_IND_PILOTS - 1)' * D / 64));
ch_21 = mean(pilots_21(:) .* exp(-1i * 2 * pi * (SC_IND_PILOTS - 1)' * D / 64));
ch_22 = mean(pilots_22(:) .* exp(-1i * 2 * pi * (SC_IND_PILOTS - 1)' * D / 64));

H_hat = [[ch_11, ch_12]; [ch_21, ch_22]];
% H_hat = [[1.25-1.55i, 0.96+1.34i], [-1.76+1.42i, -0.99-1.10i]];

syms_eq = Y_data;

ch_12 = 0;
ch_22 = 0;

phase_err = reshape(repmat(exp(1i * 2 * pi * (0:1:63)' .* D / 64), [1, 44]), 1, 64, 44);

syms_eq_pc(1, :, :) = phase_err .* syms_eq(1, :, :);
syms_eq_pc(2, :, :) = phase_err .* syms_eq(2, :, :);

syms_eq_pc(1, :, :) = (conj(ch_11 + ch_12) / ((abs(ch_11) ^ 2) + (abs(ch_12) ^ 2))) .* syms_eq_pc(1, :, :);
syms_eq_pc(2, :, :) = (conj(ch_21 + ch_22) / ((abs(ch_21) ^ 2) + (abs(ch_22) ^ 2))) .* syms_eq_pc(2, :, :);

% Take only data SCs
syms_eq_pc = syms_eq_pc(:,SC_IND_DATA,:);

% Reshape the 3-D matrix to 2-D:
syms_eq_pc = reshape(syms_eq_pc, N_UE, [] );
syms_eq_pc = syms_eq_pc.';


%% Demodulate
rx_data = demod_sym(syms_eq_pc ,MOD_ORDER);

%% Plot results

cf = 0;
fst_clr = [0, 0.4470, 0.7410];
sec_clr = [0.8500, 0.3250, 0.0980];

% Rx signal
cf = cf + 1;
figure(cf); clf;
for sp = 1:N_BS_NODE
    subplot(N_BS_NODE,2,2*(sp -1) + 1 );
    plot(real(rx_vec_iris(sp,:)));
    axis([0 length(rx_vec_iris(sp,:)) -TX_SCALE TX_SCALE])
    grid on;
    title(sprintf('BS antenna %d Rx Waveform (I)', sp));

    subplot(N_BS_NODE,2,2*sp);
    plot(imag(rx_vec_iris(sp,:)), 'color' , sec_clr);
    axis([0 length(rx_vec_iris(sp,:)) -TX_SCALE TX_SCALE]);
    grid on;
    title(sprintf('BS antenna %d Rx Waveform (Q)', sp));
end 

if PLOT
%Tx signal
cf = cf + 1;
figure(cf); clf;
for sp=1:N_UE
    subplot(N_UE,2,2*(sp -1) + 1);
    plot(real(tx_vecs_iris(:,sp)));
    axis([0 length(tx_vecs_iris(:,sp)) -TX_SCALE TX_SCALE])
    grid on;
    title(sprintf('UE %d Tx Waveform (I)', sp));

    subplot(N_UE,2,2*sp);
    plot(imag(tx_vecs_iris(:,sp)), 'color' , sec_clr);
    axis([0 length(tx_vecs_iris(:,sp)) -TX_SCALE TX_SCALE])
    grid on;
    title(sprintf('UE %d Tx Waveform (Q)',sp));
end

%% Rx LTS correlation
cf = cf+ 1;
figure(cf); clf;
for sp = 1:N_BS_NODE
    subplot(N_BS_NODE,1,sp);
    plot(lts_corr(sp,:)) 
    grid on;
    xlabel('Samples');
    y_label = sprintf('Anetnna %d',sp);
    ylabel(y_label);
    myAxis = axis();
    axis([1, 1000, myAxis(3), myAxis(4)])
end
tb = annotation('textbox', [0 0.87 1 0.1], ...
    'String', 'LTS Correlations', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
tb.FontWeight = 'bold';

%% Constellations
cf = cf+ 1;
figure(cf); clf;
if N_BS_NODE >=4
    sp_rows = ceil(N_BS_NODE/4)+1;
else
    sp_rows = ceil(N_BS_NODE/2)+1;
end
sp_cols = ceil(N_BS_NODE/(sp_rows -1));

for sp=1:N_UE
    subplot(sp_rows,sp_cols, sp);
    plot(syms_eq_pc(:,sp),'o','MarkerSize',1, 'color', sec_clr);
    axis square; axis(1.5*[-1 1 -1 1]);
    grid on;
    hold on;
    plot(tx_syms(:, 1),'*', 'MarkerSize',10, 'LineWidth',2, 'color', fst_clr);
    title(sprintf('Equalized Uplink Tx (blue) and \n Rx (red) symbols for stream %d', sp));
    % legend({'Rx','Tx'},'Location','EastOutside', 'fontsize', 12);
end

for sp=1:N_BS_NODE
    subplot(sp_rows,sp_cols, sp_cols+sp);
    plot(squeeze(Y_data(sp,:,:)),'co','MarkerSize',1);
    axis square; axis(max(max(max( abs( Y_data)) ))*[-1 1 -1 1]);
    title(sprintf('Unequalized received symbols \n at BS ant. %d', sp));
    grid on;
    hold on;
end

% %% Channel Estimates
% cf = cf + 1;
% cond_clr = [0.8500, 0.250, 0.1980];
% 
% bw_span = (20/N_SC) * (-(N_SC/2):(N_SC/2 - 1)).';
% 
% figure(cf); clf;
% sp = 0;
% for ibs = 1:N_BS_NODE
%     for iue = 1:N_UE
%         sp = sp+1;
%         subplot(N_BS_NODE,N_UE,sp); 
%         bar(bw_span, fftshift(abs( squeeze(H_hat(ibs, iue, : ) ) ) ),1,'LineWidth', 1);
%         axis([min(bw_span) max(bw_span) 0 1.1*max(abs( squeeze(H_hat(ibs, iue, :) ) ) )])
%         grid on;
%         title(sprintf('UE %d -> BS ant. %d Channel Est. (Mag.)', iue, ibs))
%         xlabel('Baseband Frequency (MHz)')
%     end
% end
% 
% subplot(N_BS_NODE+1,1,N_BS_NODE+1);
% bh = bar(bw_span, fftshift(channel_condition_db) ,1, 'LineWidth', 1);
% set(bh,'FaceColor',cond_clr);
% axis([min(bw_span) max(bw_span) 0 max(channel_condition_db)+1])
% grid on;
% title('Channel Condition (dB)')
% xlabel('Baseband Frequency (MHz)')

end

%% EVM & SNR
% tx_data = [tx_data, tx_data];
rx_data = rx_data(:, 1);
sym_errs = sum(tx_data(:) ~= rx_data(:));
bit_errs = length(find(dec2bin(bitxor(tx_data(:), rx_data(:)),8) == '1'));
evm_mat = double.empty();
aevms = zeros(N_UE,1);
snr_mat = zeros(N_UE,1);

cf = cf + 1;
figure(cf); clf;
for sp = 1:N_UE
    tx_vec = tx_syms_mat(:,:,1);
    evm_mat(:,sp)  = abs(tx_vec(:) - syms_eq_pc(:,sp) ).^2;
    aevms(sp) = mean(evm_mat(:,sp));
    snr_mat(sp) = -10*log10(aevms (sp));
    
    subplot(2,N_UE,sp)
    plot(100*evm_mat(:,sp),'o','MarkerSize',1)
    axis tight
    hold on
    plot([1 length(evm_mat(:,sp) )], 100*[aevms(sp), aevms(sp)],'color', sec_clr,'LineWidth',4)
    hold off
    xlabel('Data Symbol Index')
    ylabel('EVM (%)');
    legend('Per-Symbol EVM','Average EVM','Location','NorthWest');
    
    h = text(round(.05*length(evm_mat(:,sp))), 100*aevms(sp), sprintf('Effective SINR: %.1f dB', snr_mat(sp)));
    set(h,'Color',[1 0 0])
    set(h,'FontWeight','bold')
    set(h,'FontSize',10)
    set(h,'EdgeColor',[1 0 0])
    set(h,'BackgroundColor',[1 1 1])
    
    title(sprintf('Stream from UE %d', sp));
    grid on
    
end


for sp=1:N_UE
   subplot(2,N_UE,N_UE+sp);
   imagesc(1:N_OFDM_SYM, (SC_IND_DATA - N_SC/2), 100*fftshift( reshape(evm_mat(:,sp), [], N_OFDM_SYM), 1));
    
    hold on;
    h = line([1,N_OFDM_SYM],[0,0]);
    set(h,'Color',[1 0 0]);
    set(h,'LineStyle',':');
    hold off;
    grid on
    xlabel('OFDM Symbol Index');
    ylabel('Subcarrier Index');
    title(sprintf('Stream from UE %d', sp));
    h = colorbar;
    set(get(h,'title'),'string','EVM (%)'); 
end

fprintf('\n MIMO Results:\n');
fprintf('Num Bytes:   %d\n', N_UE*N_DATA_SYMS * log2(MOD_ORDER) / 8);
fprintf('Sym Errors:  %d (of %d total symbols)\n', sym_errs, N_UE * N_DATA_SYMS);
fprintf('Bit Errors:  %d (of %d total bits)\n', bit_errs, N_UE*N_DATA_SYMS * log2(MOD_ORDER));

fprintf("SNRs (dB): \n");
fprintf("%.2f\t", snr_mat);
fprintf("\n");
