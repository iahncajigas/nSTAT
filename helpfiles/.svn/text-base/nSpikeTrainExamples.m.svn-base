%% Test the nspikeTrain Class

%% Example 1: Using the nspikeTrain Class
% Lets create some pseudo data and plot it.

spikeTimes = sort(rand(1,100))*1;
spikeTimes = unique(round(spikeTimes*10000)./10000); %round off;
nst=nspikeTrain(spikeTimes,'n1',.001,0,1);
figure; nst.plot;

%%
% We can now change the signal representation of the nspikeTrain and see
% what effects it has.
%
% 100ms bins from 0 to 10 sec. Actual SignalObj representation of the
% nspikeTrain is not changed because are using getSigRep
figure; nst.resample(1/.1);
nst.getSigRep.plot; 

%%
% 10ms bins from 0 to 10 sec. Actually changing the representation of the
% signal.
figure; nst.resample(1/.01);
nst.getSigRep.plot; 

%%
% Get the largest binsize that still maintains a binary signal
% representation
figure; nst.resample(1/nst.getMaxBinSizeBinary);
nst.getSigRep.plot;