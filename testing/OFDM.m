clc;
clear all;

x = randn(256, 10);
y = fftshift(x);
y = ifft(y);
y = [y(end-64+1:end, :); y];

x_ = y(65:end, :);
x_ = fft(x_);
x_ = fftshift(x_);

10;