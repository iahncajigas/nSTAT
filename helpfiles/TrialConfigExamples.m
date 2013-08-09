%% TrialConfig Examples
% tcObj=TrialConfig(covMask,sampleRate, history,minTime,maxTime)
tc1 = TrialConfig({'Force','f_x'},2000,[.1 .2],-1,2);
tc2 = TrialConfig({'Position','x'},2000,[.1 .2],-1,2);
tcc = ConfigColl({tc1,tc2});