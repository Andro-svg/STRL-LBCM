close all;
clear all;
clc
clear

format long
addpath('utils/');
addpath('metric_utils\');

frame = 16;
lambdaL = 10;
mu = 0.05; 

readPath = '.\data';
savePath = '.\result';  

if ~exist(savePath)
    mkdir(savePath);
end

tuneopts.temporal_step = frame;
tuneopts.lambdaL = lambdaL;
tuneopts.mu = mu;
target_detection(char(readPath), savePath, tuneopts);

