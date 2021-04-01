
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Author(s): Yashwanth R - yashwanthr@iisc.ac.in
%
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
WRITE_PNG_FILES         = 1;           % Enable writing plots to PNG
SIM_MOD                 = 0;
% DEBUG                   = 1;
PLOT                    = 1;

keepTXRX_Running = 1;

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
    
    bs_ids = ["RF3E000064", "RF3E000037", "RF3E000012"];
    numBSAnts = 16;

    ue_ids= ["RF3E000044"];
    numUEAnts = 16;

    N_BS_NODE               = length(bs_ids);           % Number of nodes/antennas at the BS
    N_UE                    = length(ue_ids);           % Number of UE nodes

    numAntsPerBSNode = numBSAnts / N_BS_NODE;
    numAntsPerUENode = numUEAnts / N_UE;

end
fprintf("5G Testbed@ IISc: IRIS MIMO Setup\n");
fprintf("Transmission type: %s \n",chan_type);

while (1)
    
    %% Downlink Baseband Signal Generation
    N = 256;
    CP = N/4;
    
    codeRate = 1/2;
    modOrder = 4;
    
    PILOTS_SYMS_IND = [1;6;11];
    DATA_SYMS_IND = [2;3;4;5;7;8;9;10];

    nonZeroSubcarriers = setdiff(1:N, []);
    
    numBits = (1/codeRate) * (1/log2(modOrder)) * length(DATA_SYMS_IND) * N;

    numSyms = floor(4096 / (N + CP));
    
    tx_data_buff = zeros(N + CP, numSyms - 1);
    
    data = randi([0, 1], 
    
    
    
    %% Burn Data onto IRIS BS Nodes
    
    
    %% Retreive Data from IRIS UE Antenna Buffers

    
    %% Downlink Baseband Signal Processing at IRIS UE Nodes
    
    
    %% Evaluation 
    
    
    if (keepTXRX_Running == 0)
        break;
    end
end


