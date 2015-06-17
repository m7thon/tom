function [data, nO, nU] = read_sequence(filename)
% reads a Sequence in tom format from the given filename
addpath('jsonlab');
data = loadjson(filename);
nU = data.nU;
nO = data.nO;
data = data.data;
