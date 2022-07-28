clc
clear all
close all

NSample = 5000;
N = 100;
normVec = zeros(NSample,1);
sig_max = 0.01;
sig_min = 0.001;
for i = 1:NSample
%     sig_i = sig_min * (sig_max-sig_min)*rand(1,1);
%     r = normrnd(zeros(N,1),sig_i*ones(N,1));
    r = normrnd(zeros(N,1),0.001*ones(N,1));
    normVec(i) = norm(r);
end

plot(normVec)