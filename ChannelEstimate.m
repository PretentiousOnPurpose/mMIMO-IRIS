function Hhat = ChannelEstimate(rx_data_buff, t1, PILOTS_SYMS_IND, DATA_SYMS_IND, nonZeroSubcarriers)

% Extract Tansmitted Pilot Symbol
tx_pilot_buff = t1(nonZeroSubcarriers, PILOTS_SYMS_IND);

% Extract Received Pilot Symbol
rx_pilot_buff = rx_data_buff(nonZeroSubcarriers, PILOTS_SYMS_IND);

% Least Square Channel Estimation
HLS = rx_pilot_buff./tx_pilot_buff;

Hhat = zeros(size(rx_data_buff));
Hr = zeros(numel(nonZeroSubcarriers), 11);
Hi = Hr;

% Interpolate
for k = 1:numel(nonZeroSubcarriers)
    Hr(k, :) = interp1(PILOTS_SYMS_IND, real(HLS(k,:)), union(PILOTS_SYMS_IND,DATA_SYMS_IND))';
    Hi(k, :) = interp1(PILOTS_SYMS_IND, imag(HLS(k,:)), union(PILOTS_SYMS_IND,DATA_SYMS_IND))';
end

H = complex(Hr, Hi);

Hhat(nonZeroSubcarriers,:) = H;