clc
clear all
close all


LoadingPath = '/Users/nzcc-ghani/Documents/MATLAB/SensorConnectivity/DATA-Sensor-GY';


s =10; %Number of subjects
DataName = 'DATA-Sensor-GY';
AvgEasyData = zeros(17,17);
AvgMediumData = zeros(17,17);

AvgHardData = zeros(17,17);


for i = 1:s
 %% Setting up paths and PLV data saved from brainstorm  
        EasyPath =[LoadingPath '/Participant ' num2str(i) '/Easy']; 
        MediumPath =[LoadingPath '/Participant ' num2str(i) '/Medium']; 
        HardPath =[LoadingPath '/Participant ' num2str(i) '/Hard']; 

% Define the path and half of the file name
%% Load data
EasyData = load([EasyPath '/DMN.mat']);
EasyData = EasyData.DMN;

MediumData = load([MediumPath '/DMN.mat']);
MediumData = MediumData.DMN;

HardData = load([HardPath '/DMN.mat']);
HardData = HardData.DMN;
%% Calculating average

AvgEasyData = AvgEasyData + EasyData;
AvgMediumData = AvgMediumData + MediumData;
AvgHardData = AvgHardData + HardData;
end

AvgEasyData = AvgEasyData/s;
AvgMediumData = AvgMediumData/s;
AvgHardData = AvgHardData/s;

save ConnMatrices AvgEasyData AvgMediumData AvgHardData


