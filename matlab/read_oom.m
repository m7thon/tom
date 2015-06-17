function [sig, tau, w0] = read_oom(filename)
% reads an Oom in tom format from the given filename
% tau will be an array of size dim x dim x nO for OOMs or of size
% dim x dim x nO x nU for IO-OOMs
addpath('jsonlab');
data = loadjson(filename);
sig = data.sig;
w0 = data.w0;
if data.nU == 0
    tau = reshape(transpose(data.tau), data.dim, data.dim, data.nO);
    for o = 1:data.nO
        tau(:,:,o) = transpose(tau(:,:,o));
    end
else
    tau_intermediary = reshape(transpose(data.tau), data.dim, data.dim, data.nU, data.nO);
    tau = zeros(data.dim, data.dim, data.nO, data.nU);
    for a = 1:data.nU
        for o = 1:data.nO
            tau(:,:,o, a) = transpose(tau_intermediary(:,:,a,o));
        end
    end
end