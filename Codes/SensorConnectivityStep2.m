clc
clear all
close all


LoadingPath = '/Users/nzcc-ghani/Documents/MATLAB/SensorConnectivity/';


s =26; %Number of subjects
DataName = 'DATA-Sensor-GZ';


     network = 'DMN'; 
     frequency = 'Alpha'; 

    ln = 17;

restoredefaultpath
addpath /Users/nzcc-ghani/Documents/FieldTrip
ft_defaults
    
    %% Baseline vs. left stimulation- data struct
    
       
        EasyVsMedium(1).label = {['Time segment ' num2str(1)]};
        EasyVsMedium(1).freq = 1:ln;
        EasyVsMedium(1).time = 1:ln;
        EasyVsMedium(1).dimord = 'rpt_chan_freq_time';
        EasyVsMedium(1).powspctrm = zeros(s, 1, ln, ln);
        %Post
        EasyVsMedium(2).label = {['Time segment ' num2str(1)]};
        EasyVsMedium(2).freq = 1:ln;
        EasyVsMedium(2).time = 1:ln;
        EasyVsMedium(2).dimord = 'rpt_chan_freq_time';
        EasyVsMedium(2).powspctrm = zeros(s, 1, ln, ln);
%% Filling EasyVsMedium data struct    
     for i= 1:s
        
       EasyPath =[LoadingPath DataName '/Participant '  num2str(i) '/Easy']; 
       MediumPath =[LoadingPath DataName '/Participant '  num2str(i) '/Medium']; 
       HardPath =[LoadingPath DataName '/Participant '  num2str(i) '/Hard']; 
         
          
        DataE0(i)=load([EasyPath '/DMN.mat'],network);
        EasyVsMedium(2).powspctrm(i, 1, :, :) =  DataE0(i).(network);  % this should yield a 1 x 12 x 5 matrix

        DataM0(i)=load([MediumPath '/DMN.mat'], network);
        EasyVsMedium(1).powspctrm(i, 1, :, :) =  DataM0(i).(network);
    end

    
    %% Filling Easy vs Hard data struct     
        EasyVsHard(1).label = {['Time segment ' num2str(1)]};
        EasyVsHard(1).freq = 1:ln;
        EasyVsHard(1).time = 1:ln;
        EasyVsHard(1).dimord = 'rpt_chan_freq_time';
        EasyVsHard(1).powspctrm = zeros(s, 1, ln, ln);
        
        EasyVsHard(2).label = {['Time segment ' num2str(1)]};
        EasyVsHard(2).freq = 1:ln;
        EasyVsHard(2).time = 1:ln;
        EasyVsHard(2).dimord = 'rpt_chan_freq_time';
        EasyVsHard(2).powspctrm = zeros(s, 1, ln, ln);
     
        for i= 1:s
        
        EasyPath =[LoadingPath DataName '/Participant '  num2str(i) '/Easy']; 
        MediumPath =[LoadingPath DataName '/Participant '  num2str(i) '/Medium']; 
        HardPath =[LoadingPath DataName '/Participant '  num2str(i) '/Hard']; 
         
          
        DataE1(i)=load([EasyPath '/DMN.mat'], network);
        EasyVsHard(2).powspctrm(i, 1, :, :) =  DataE1(i).(network);  % this should yield a 1 x 12 x 5 matrix

        DataH1(i)=load([HardPath '/DMN.mat'], network);
        EasyVsHard(1).powspctrm(i, 1, :, :) =  DataH1(i).(network);
         end
%  
%        
    %% Left vs. right stimulation- data struct
% 
        MediumVsHard(1).label = {['Time segment ' num2str(1)]};
        MediumVsHard(1).freq = 1:ln;
        MediumVsHard(1).time = 1:ln;
        MediumVsHard(1).dimord = 'rpt_chan_freq_time';
        MediumVsHard(1).powspctrm = zeros(s, 1, ln, ln);
        
        MediumVsHard(2).label = {['Time segment ' num2str(1)]};
        MediumVsHard(2).freq = 1:ln;
        MediumVsHard(2).time = 1:ln;
        MediumVsHard(2).dimord = 'rpt_chan_freq_time';
        MediumVsHard(2).powspctrm = zeros(s, 1, ln, ln);
     
        for i= 1:s
        
        EasyPath =[LoadingPath DataName '/Participant '  num2str(i) '/Easy']; 
        MediumPath =[LoadingPath DataName '/Participant '  num2str(i) '/Medium']; 
        HardPath =[LoadingPath DataName '/Participant '  num2str(i) '/Hard']; 
        
          
        DataM2(i)=load([MediumPath '/DMN.mat'], network);
        MediumVsHard(2).powspctrm(i, 1, :, :) =  DataM2(i).(network);  % this should yield a 1 x 12 x 5 matrix

        DataH2(i)=load([HardPath '/DMN.mat'], network);
        MediumVsHard(1).powspctrm(i, 1, :, :) =  DataH2(i).(network);
         end

    
    %% Config for freqstatistics calculations
    cfg = [];
    cfg.method            = 'montecarlo';           % using the Monte Carlo Method to calculate the significance probability
    cfg.statistic         = 'depsamplesT';        % using the independent samples T-statistic as a measure to evaluate the effect at the sample level
    cfg.correctm          = 'cluster';
    cfg.clusteralpha      = 0.1;                   % alpha level of the sample-specific test statistic that will be used for thresholding
    cfg.clustertail       = 0;
    cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
    cfg.tail              = 0;                      % -1, 1 or 0; one-sided or two-sided test
    cfg.correcttail       = 'prob';                 % the two-sided test imPLIes that we do non-parametric two tests
    cfg.alpha             = 0.05;                   % alpha level of the permutation test
    cfg.numrandomization  = 5000;                   % number of draws from the permutation distribution
    cfg.neighbours        = [];                    % number of draws from the permutation distribution
    cfg.design            = [ones(1,s) 2*ones(1,s); 1:s 1:s];% design matrix
    cfg.ivar              = 1;                      % the index of the independent variable in the design matrix
    cfg.uvar              = 2;
   % cfg.neighbours        = [];
    %cfg.numrandomization= 'all'; % there are no spatial neighbours, only in time and frequency
    % run sourcestatistics using cluster based correction %
  % Use any of the data matrices


    
    %% Freqstatistics calculations
    % Baseline vs. left stimulation
   StatEasyVsMed = ft_freqstatistics(cfg, EasyVsMedium(2), EasyVsMedium(1));

    % Baseline vs. right stimulation
   StatEasyVsHard = ft_freqstatistics(cfg, EasyVsHard(2), EasyVsHard(1));
        
    % Left vs. right stimulation
   StatMedVsHard = ft_freqstatistics(cfg, MediumVsHard(2), MediumVsHard(1));
   
   %% Changing dimentions of probability and t values
t_EvH = squeeze(StatEasyVsHard.stat);
p_EvH = squeeze(StatEasyVsHard.prob);

t_EvM = squeeze(StatEasyVsMed.stat);
p_EvM = squeeze(StatEasyVsMed.prob);

t_MvH = squeeze(StatMedVsHard.stat);
p_MvH = squeeze(StatMedVsHard.prob);
   
   %% Calculate significant areas (Easy vs Medium)


        Med = EasyVsMedium(2).powspctrm(1:s,1,:,:);              % select the Medium trial
        Easy = EasyVsMedium(1).powspctrm(1:s,1,:,:);              % select the Easy trials
   
         MeanPS_Medium  = mean(Med);
         MeanPS_Easy = mean(Easy);

    MeanPS_Medium_std  = std(Med);
    MeanPS_Easy_std = std(Easy);
    Size1    = size(MeanPS_Medium);
%% Calculating size that will reshape that data after effect size calculation
BS = MeanPS_Easy - MeanPS_Medium;
Size = size(BS);
%% Reducing dimension of the data from 4D to calculate effect size
MeanPS_Medium = squeeze(reshape(MeanPS_Medium, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Easy = squeeze(reshape(MeanPS_Easy, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Medium_std = squeeze(reshape(MeanPS_Medium_std, Size1(2:end)));
MeanPS_Easy_std = squeeze(reshape(MeanPS_Easy_std, Size1(2:end)));
%% Calculating effect size
n1 = size(MeanPS_Medium,1);
n2 = size(MeanPS_Easy,1);
    pooled_sd = sqrt(((n1-1)*(MeanPS_Medium_std^2) + (n2-1)*(MeanPS_Easy_std^2))/(n1+n2-2));

% Calculate Cohen's d
EffectSizeEasyVsMedium1 = (MeanPS_Medium - MeanPS_Easy) / pooled_sd;
%% reshape data and assign it to stat structure for plotting
EffectSizeEasyVsMedium = reshape(EffectSizeEasyVsMedium1,Size(2:end));
StatEasyVsMed.effect = EffectSizeEasyVsMedium;
%   

% StatEasyVsMed.prob(isnan(StatEasyVsMed.prob)) = 1;

      %% Calculate significant areas (Easy vs Hard)


        Hard = EasyVsHard(2).powspctrm(1:s,1,:,:);              % select the Medium trial
        Easy = EasyVsHard(1).powspctrm(1:s,1,:,:);              % select the Easy trials
   
         MeanPS_Hard  = mean(Hard);
         MeanPS_Easy = mean(Easy);

    MeanPS_Hard_std  = std(Hard);
    MeanPS_Easy_std = std(Easy);
    Size1    = size(MeanPS_Hard);
%% Calculating size that will reshape that data after effect size calculation
BS = MeanPS_Easy - MeanPS_Hard;
Size = size(BS);
%% Reducing dimension of the data from 4D to calculate effect size
MeanPS_Hard = squeeze(reshape(MeanPS_Hard, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Easy = squeeze(reshape(MeanPS_Easy, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Hard_std = squeeze(reshape(MeanPS_Hard_std, Size1(2:end)));
MeanPS_Easy_std = squeeze(reshape(MeanPS_Easy_std, Size1(2:end)));
%% Calculating effect size
n1 = size(MeanPS_Hard,1);
n2 = size(MeanPS_Easy,1);
    pooled_sd = sqrt(((n1-1)*(MeanPS_Hard_std^2) + (n2-1)*(MeanPS_Easy_std^2))/(n1+n2-2));

% Calculate Cohen's d
EffectSizeEasyVsHard1 = (MeanPS_Easy - MeanPS_Hard) / pooled_sd;
%% reshape data and assign it to stat structure for plotting
EffectSizeEasyVsHard = reshape(EffectSizeEasyVsHard1,Size(2:end));
StatEasyVsHard.effect = EffectSizeEasyVsHard;
%   

StatEasyVsHard.prob(isnan(StatEasyVsHard.prob)) = 1;

  %% Calculate significant areas ( Medium vs Hard)


        Med = MediumVsHard(2).powspctrm(1:s,1,:,:);              % select the Medium trial
        Hard = MediumVsHard(1).powspctrm(1:s,1,:,:);              % select the Easy trials
   
         MeanPS_Medium  = mean(Med);
         MeanPS_Hard = mean(Hard);

    MeanPS_Medium_std  = std(Med);
    MeanPS_Hard_std = std(Hard);
    Size1    = size(MeanPS_Medium);
%% Calculating size that will reshape that data after effect size calculation
BS = MeanPS_Hard - MeanPS_Medium;
Size = size(BS);
%% Reducing dimension of the data from 4D to calculate effect size
MeanPS_Medium = squeeze(reshape(MeanPS_Medium, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Hard = squeeze(reshape(MeanPS_Hard, Size1(2:end))); % we need to "squeeze" out one of the dimensions, i.e. make it 3-D rather than 4-D
MeanPS_Medium_std = squeeze(reshape(MeanPS_Medium_std, Size1(2:end)));
MeanPS_Hard_std = squeeze(reshape(MeanPS_Hard_std, Size1(2:end)));
%% Calculating effect size
n1 = size(MeanPS_Medium,1);
n2 = size(MeanPS_Hard,1);
    pooled_sd = sqrt(((n1-1)*(MeanPS_Medium_std^2) + (n2-1)*(MeanPS_Hard_std^2))/(n1+n2-2));

% Calculate Cohen's d
EffectSizeMediumVsHard1 = (MeanPS_Hard - MeanPS_Medium) / pooled_sd;
% EffectSizeMediumVsHard1(EffectSizeMediumVsHard1>0) = 0;

%% reshape data and assign it to stat structure for plotting
EffectSizeMediumVsHard = reshape(EffectSizeMediumVsHard1,Size(2:end));

% EffectSizeMediumVsHard1(EffectSizeMediumVsHard1>0) = 0;

StatMedVsHard.effect = EffectSizeMediumVsHard;
%   

StatMedVsHard.prob(isnan(StatMedVsHard.prob)) = 1;


%%
    %% Config for stat plot
    cfg.renderer      = 'openGL';     % painters does not support opacity, openGL does
    cfg.colorbar      = 'yes';
    cfg.parameter     = 'effect'; % display the statistical value, i.e. the t-score
    cfg.maskparameter = 'mask';       % use significance to mask the power
    cfg.maskalpha     = 0;          % make non-significant regions 30% visible
    cfg.zlim = [-20 0];

   cfg.interactive    =  'yes';
%    cfg.directionality = 'inflow';
    
    figure()
    %set(gcf,'name',['Statistics' network ': Within, Control'],'NumberTitle','on')
    ft_singleplotTFR(cfg, StatEasyVsMed);
    title(['significant power changes (p<0.05)' frequency])
    ylabel('Brain region')
    xlabel('Brain region')
    
    figure()
    %set(gcf,'name',['Statistics ' network ': Within, Chiro'],'NumberTitle','on')
    ft_singleplotTFR(cfg, StatEasyVsHard);
    title(['significant power changes (p<0.05)' frequency])
    ylabel('Brain region')
    xlabel('Brain region')
    
     figure()
%     %set(gcf,'name',['Statistics ' network ': Within, Chiro'],'NumberTitle','on')
     ft_singleplotTFR(cfg, StatMedVsHard);
     title(['significant power changes (p<0.05)' frequency])
     ylabel('Brain region')
     xlabel('Brain region')




