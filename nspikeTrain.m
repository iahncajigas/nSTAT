classdef nspikeTrain < handle
% NSPIKETRAIN A neural spike train object consists of a sequence of
% spikeTimes. The spike train can be represented as a SignalObj with a
% particular sampling rate. The 1/sampleRate (sampling period) is larger
% that the difference between any two spike times, the neural spike train
% will no longer have a binary representation.
%
% Usage:
% nst=nspikeTrain(spikeTimes,name,binwidth,minTime,maxTime, varargin)
%      spikeTimes: row or column vector of spike times. 
%
%      OPTIONAL INPUTS:
%      name:       name of neuron data recorded from. Default='';
%
%      binwidth:   binwidth to be used for SignalObj representation of
%                  spikeTimes. Default: 0.01 sec/bin
%
%      minTime:    Default is min(spikeTimes)
%
%      maxTime:    Default is max(spikeTimes)
%
%      varargin: xlabelval, xunits, yunits,dataLabels in that order can be
%      passed to the SignalObj constructor for the signal representation of
%      the neural spike train.
%
% <a href="matlab: methods('nspikeTrain')">methods</a>
% <a href="matlab:web('nSpikeTrainExamples.html', '-helpbrowser')">nSpikeTrain Examples</a> 
%
% see also <a href="matlab:help('Covariate')">Covariate</a>, <a 
% href="matlab:help('SignalObj')">SignalObj</a>
% 
% Reference page in Help browser
% <a href="matlab: doc('nspikeTrain')">doc nspikeTrain</a>


%
% nSTAT v1 Copyright (C) 2012 Masschusetts Institute of Technology
% Cajigas, I, Malik, WQ, Brown, EN
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License 
% along with this program; if not, write to the Free Software Foundation, 
% Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

    properties (SetAccess = private)
        name        % name of the nspikeTrain
        spikeTimes  % collection of times at which spikes occured
        
        sigRep      % SignalObj representation of nspikeTrain 
        sampleRate  % sampleRate for the sigRep
        maxTime     % maximum time of interest or time of last spike
        minTime     % minimum time of interest or time of first spike
        
        %%TODO add listener to each spike train so that consistency is
        %%guaranteed of objects are modified.
        isSigRepBin % Boolean indicating 1 or 0 spikes occur per bin
        xlabelval
        xunits
        yunits
        dataLabels
        MER  %signal of the micro-electrode recordings when available 
        avgFiringRate
        %Bursting related parameters%
        B    %burstiness parameter
        An   %burtiness parameters without finite sample effect - Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite event sequences. Physical Review E, 94(3). http://doi.org/10.1103/physreve.94.032311
        burstTimes
        burstRate
        burstDuration
        burstSig
        burstIndex % Hutchinson, et al. 1997 Effects of Apomorphine on GP Neurons in parkinsian Patients
        numBursts
        numSpikesPerBurst
        avgSpikesPerBurst
        stdSpikesPerBurst
        Lstatistic  % Goldberg et. al 2002. Enhanced Synchrony among Primary Motor Cortex Neurons in the
                    % 1-Methyl-4-Phenyl-1,2,3,6-Tetrahydropyridine Primate Model of Parkinson’s Disease
    end
    
    methods
        function nst=nspikeTrain(spikeTimes,name,binwidth,minTime,maxTime, xlabelval, xunits, yunits,dataLabels,makePlots)
            %constructor
            if(nargin<10|| isempty(makePlots))
                makePlots=0;
            end
            if(nargin<9 || isempty(dataLabels))
                dataLabels = '';
            end
            if(nargin<8 || isempty(yunits))
                yunits='';
            end
            if(nargin<7 || isempty(xunits))
                xunits='s';
            end
            if(nargin<6 || isempty(xlabelval))
                xlabelval='time';
            end
            if(nargin<5)
                maxTime = max(spikeTimes);
                if(isempty(maxTime))
                    maxTime=0;
                end
            end
            if(nargin<4)
                minTime= min(spikeTimes);
                if(isempty(minTime))
                    minTime=0;
                end
            end
            if(nargin<3)
                binwidth=.001;
            end
            if(nargin<2)
                name='';
            end
            if(nargin<1)
                error('nspikeTrain requires a spikeTimes array as input to create an object');
            end
            [l,w]=size(spikeTimes);
            if(l>w)
                nst.spikeTimes=spikeTimes';      
            else
                nst.spikeTimes=spikeTimes;
            end
            
            nst.name=name;
            nst.sampleRate = 1/binwidth;
            nst.minTime = minTime;
            nst.maxTime = maxTime;
            nst.xlabelval = xlabelval;
            nst.xunits = xunits;
            nst.yunits = yunits;
            nst.dataLabels = dataLabels;
%             nst.setSigRep(binwidth, minTime, maxTime,varargin{:});

            nst.sigRep = [];
            nst.isSigRepBin=[];

            nst.computeStatistics(makePlots);            
        end
        function Lstat = getLStatistic(nstObj)
           mISIs = mean(nstObj.getISIs); 
           Pt = nstObj.getSigRep(mISIs);
           Lstat = length(unique(Pt.data));
         
        end
%         function shift(nstObj,deltaT)
%            nstObj.spikeTimes = nstObj.spikeTimes + deltaT;
%            nstObj.setMinTime(nstObj.minTime+deltaT);
%            nstObj.setMaxTime(nstObj.maxTime+deltaT);
%            nstOb.sigRep = [];
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setMER(nstObj,MERSig)
           if(isa(MERSig,'SignalObj'))
                nstObj.MER = MERSig;
           end
        end
           
        function setName(nstObj,name)
           % setName(nstObj,name) 
           % set the name after construction
           nstObj.name = name;
%            if(~isempty(nstObj.sigRep))
%                 nstObj.sigRep.setName(name);
%            end
        end
        
        function computeStatistics(nstObj,makePlots)
            if(nargin<2 || isempty(makePlots))
                makePlots = 0;
            end
            ISI=nstObj.getISIs;
            spikeTimes =nstObj.spikeTimes;    
            nstObj.avgFiringRate = length(spikeTimes)/(nstObj.maxTime - nstObj.minTime);
            
            
            %% Compute Burst Parameters
            % Hutchinson, et al. 1997
            if(~isempty(ISI))
                nstObj.burstIndex=1/mode(ISI)/nstObj.avgFiringRate;
            
                % Cells havingburst index values of approximately 1 to 1.5 had regular firing 
                % intervals (ie, Border cells), cells in the range of 1.6 to 9.9 were “irregular” firing cells,
                % and those with values greater that 10 were termed bursting cells.
            
            
                % Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite 
                % event sequences. Physical Review E, 94(3), 032311. 
                % http://doi.org/10.1103/PhysRevE.94.032311
                sigma=std(ISI);
                mu = mean(ISI);
                r= sigma/mu;
                nstObj.B = (r-1)/(r+1);  %burstiness index for infinite sequence
                %B has the value of ?1 for regular time series as ? = 0, and 0 for Poissonian or random time series as ? = µ. Finally, the value of B approaches 1 for extremely bursty time series as ? ? ? for finite µ.
                n=length(spikeTimes);
                nstObj.An=(sqrt(n+2)*r-sqrt(n))./((sqrt(n+2)-2)*r+sqrt(n)); %corrected burstiness index for finite sequence

               % Chen, L., Deng, Y., Luo, W., Wang, Z., & Zeng, S. (2009). Detection of bursts in 
               % neuronal spike trains by the mean inter-spike interval method. Progress in Natural 
               % Science, 19(2), 229. http://doi.org/10.1016/j.pnsc.2008.05.027
                Ln = ISI(ISI<mu);
                ML = mean(Ln);
                if(~isempty(nstObj.spikeTimes))
                    t=spikeTimes;
    %                 t=[nstObj.spikeTimes(1); t];
                else
                    t=[];
                end
                burstISI = double(ISI<ML);
    %             B=[1 1];
    %             A=1;
    %             if(length(burstISI<4))

                    y=(burstISI(1:end)+[burstISI(2:end);0])>1;

    %             else
    % %                 y=filtfilt(B,A,burstISI)>1;
    %             end
                diffSig = [0;diff(y)];
                tdiff = (t(1:end-1)+t(2:end))/2;




                burstStart = find(diffSig==1);
                burstEnd=find(diffSig==-1)+1;
                if(isempty(burstStart))
                    burstEnd=[];
                    nstObj.burstDuration = [];
                    nstObj.burstSig = [];
                    nstObj.numSpikesPerBurst =[];
                    nstObj.numBursts=[];
                    nstObj.burstRate=[];

                end

                if(length(burstEnd)>length(burstStart))
                    burstStart = [find(y(1:burstEnd(1))==1,1, 'first'); burstStart];
                end

                if(length(burstStart)>length(burstEnd))
                    burstEnd = [find(y(burstStart(end):end)==1, 1,'last'); burstEnd];
                end

                if(~isempty(burstStart))
                    if(makePlots==1)
                        close all;
                        nstObj.plot; hold on;
                        plot(tdiff, ISI,'ko');
                        plot([t(1) t(end)], ML*[1;1]); 
                        plot(tdiff, diffSig,'r--');
                        axis tight;
                        plot(t(burstStart),1.2*ones(length(burstStart)),'bo'); hold on;
                        plot(t(burstEnd),1.2*ones(length(burstEnd)),'bd');
                    end
                    burstData = zeros(length(spikeTimes),1);
                    for i=1:length(burstStart)
                        burstData(burstStart(i):burstEnd(i))=1;
                    end

                    nstObj.burstDuration = nstObj.spikeTimes(burstEnd)-nstObj.spikeTimes(burstStart);
                    nstObj.burstSig = SignalObj(nstObj.spikeTimes,burstData,'Burst Signal');
                    nstObj.burstTimes = nstObj.spikeTimes(burstStart);
                    nstObj.numBursts=length(burstStart);
                    nstObj.burstRate=nstObj.numBursts/(nstObj.maxTime-nstObj.minTime);
                    nstObj.numSpikesPerBurst =burstEnd-burstStart + 1;
                    nstObj.avgSpikesPerBurst = mean(nstObj.numSpikesPerBurst+1);
                    nstObj.stdSpikesPerBurst = std(nstObj.numSpikesPerBurst+1);
                    nstObj.Lstatistic        = nstObj.getLStatistic; %Goldberg et al 2002.
                end
            end    
            
        end
        
        function sigRep = setSigRep(nstObj,binwidth,minTime,maxTime)
            %sigRep = setSigRep(nstObj, binwidth, minTime, maxTime)
              nstObj.sigRep = nstObj.getSigRep(binwidth,minTime,maxTime);
              nstObj.isSigRepBin = nstObj.isSigRepBinary;
              nstObj.sampleRate = nstObj.sigRep.sampleRate;
%               sigRep=nstObj.sigRep;
              nstObj.setMinTime(nstObj.sigRep.minTime);
              nstObj.setMaxTime(nstObj.sigRep.maxTime);
        end        
        function setMinTime(nstObj,minTime)
           % setMinTime(sObj,nstObj)
           % sets the minimun value of the time vector for the SignalObj representation of the nspikeTrain.
           %nstObj.sigRep.setMinTime(minTime);
           nstObj.minTime=minTime;
           nstObj.computeStatistics;
%            nstObj.clearSigRep;
        end        
%         function answer = get.isSigRepBin(nstObj)
%             answer = nstObj.isSigRepBinary;
%         end
        function setMaxTime(nstObj,maxTime)
           % setMaxTime(sObj,nstObj)
           % sets the maximum value of the time vector for the SignalObj
           % representation of the nspikeTrain.
           %nstObj.sigRep.setMaxTime(maxTime);
           nstObj.maxTime=maxTime;
           nstObj.computeStatistics;
%            nstObj.clearSigRep;
        end
        function clearSigRep(nstObj)
           nstObj.sigRep=[];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nstObj = resample(nstObj,sampleRate)
            % resample(nstObj,sampleRate)
            % change the sampleRate or equivalently the binwidth of the
            % SignalObj representation of the nSpikeTrain
            % may change the nstObj.isSigRepBin value is the binwidth is
            % great than nstObj.getMaxBinSizeBinary
%             s=nstObj.sigRep;
            binwidth=1/sampleRate;
            nstObj.setSigRep(binwidth,nstObj.minTime,nstObj.maxTime);
            nstObj.sampleRate = sampleRate;
        end       

        
        function sigRep=getSigRep(nstObj,binwidth,minTime,maxTime)
            % sigRep=getSigRep(nstObj,binwidth,minTime,maxTime)
            % returns the SignalObj representation of the nspikeTrain with
            % the parameters specified.
            %
            % binwidth: defaults to 1/nstObj.sampleRate if argument missing
            % or
            %           empty.
            %
            % minTime : defaults to nstObj.minTime if argument missing or
            %           empty.
            %
            % maxTime : defaults to nstObj.maxTime if argument missing or
            %           empty.
            %if(nargin==1)
            %    sigRep = nstObj.sigRep;
            %else
                % varargin: xlabelval, xunits, yunits,dataLabels
                if((nargin<4) || isempty(maxTime))
                    maxTime=nstObj.maxTime;
                end
                if((nargin<3) || isempty(minTime))
                    minTime=nstObj.minTime;
                end
                if((nargin<2) || isempty(binwidth))
                    binwidth=1/nstObj.sampleRate; 
                    precision =2*ceil(log10(nstObj.sampleRate));
                    binwidth = roundn(binwidth,-precision); 
                end
                precision =2*ceil(log10(1/binwidth));
                binwidth = roundn(binwidth,-precision); 
                if(and(~isempty(maxTime),~isempty(minTime)))
                %                     timeVec=linspace(minTime,maxTime,ceil((1/binwidth)*abs(maxTime-minTime)/binwidth)*binwidth+1); %scaling by binwidth to avoid roundoff error
                timeVec=linspace(minTime,maxTime,(maxTime-minTime)./binwidth +1);%:binwidth:maxTime;
                windowTimes=[minTime-binwidth/2 timeVec+binwidth/2];
                else
                  timeVec = [];
                  windowTimes=[];
                end
                data=zeros(length(timeVec),1);
                  

    
                  %If we already have the right signal representation they dont
                  %waste time.
                   if(~isempty(nstObj.sigRep))
                        if((nstObj.sigRep.sampleRate==nstObj.sampleRate) && min(nstObj.sigRep.time)==minTime && max(nstObj.sigRep.time)==maxTime)
                             sigRep = nstObj.sigRep.copySignal;
                        else %create the appropriate representation

                            spikeTimes = nstObj.spikeTimes;
%                             spikeTimes = round(spikeTimes*nstObj.sampleRate*2)/(nstObj.sampleRate*2);
                            spikeTimes = roundn(spikeTimes,-precision);
%                             windowTimes    = round(windowTimes*nstObj.sampleRate*2)/(2*nstObj.sampleRate);
                            windowTimes    = roundn(windowTimes,-precision-1);
                            lwindowTimes = length(windowTimes);
                            for j=1:length(timeVec) %number of bins
                                if(j==(lwindowTimes-1))
                                    tempSpikes=spikeTimes(spikeTimes>=windowTimes(j));
                                    data(j)=sum(tempSpikes<=windowTimes(j+1));
%                                     data(j) = sum((spikeTimes>=windowTimes(j) & spikeTimes<=windowTimes(j+1)));
                                elseif(j>floor(lwindowTimes/2))
                                    tempSpikes=spikeTimes(spikeTimes>=windowTimes(j));
                                    data(j)=sum(tempSpikes<windowTimes(j+1));
%                                     data(j) = sum((spikeTimes>=windowTimes(j) & spikeTimes<windowTimes(j+1)));
                                else
                                    tempSpikes=spikeTimes(spikeTimes<windowTimes(j+1));
                                    data(j)=sum(tempSpikes>=windowTimes(j));
                                end
                            end
                        
%                             tV=repmat(timeVec,[length(spikeTimes) 1]);
%                             sT=repmat(spikeTimes',[1 length(timeVec)-1]);
%                             data= sT>=tV(:,1:end-1) & sT<tV(:,2:end);
%                             data(:,end)= sT(:,end)>=tV(:,end-1) & sT(:,end)<=tV(:,end);
% %                             data=[sum(data) 0];
                         
                            sigRep = SignalObj(timeVec, data',nstObj.name,nstObj.xlabelval, nstObj.xunits, nstObj.yunits, nstObj.dataLabels);
                            nstObj.sigRep = sigRep;
                            if(max(sigRep.data)>1)
                                nstObj.isSigRepBin=0;
                            else
                                nstObj.isSigRepBin=1;
                            end
%                             nstObj.isSigRepBin=nstObj.isSigRepBinary;
                        end 
                   else
                        %rounding avoids comparison errors due to
                        %differences in non-significant digits
                        spikeTimes = nstObj.spikeTimes;
%                             spikeTimes = round(spikeTimes*nstObj.sampleRate*2)/(nstObj.sampleRate*2);
                        spikeTimes = roundn(spikeTimes,-precision);
%                             windowTimes    = round(windowTimes*nstObj.sampleRate*2)/(2*nstObj.sampleRate);
                        windowTimes    = roundn(windowTimes,-precision-1);
                        lwindowTimes = length(windowTimes);
                        %                         ltimeVec = length(timeVec);
                        for j=1:length(timeVec) %number of bins
                            if(j==(lwindowTimes)-1)
                                tempSpikes=spikeTimes(spikeTimes>=windowTimes(j));
                                data(j)=sum(tempSpikes<=windowTimes(j+1));
                        %                                     data(j) = sum((spikeTimes>=windowTimes(j) & spikeTimes<=windowTimes(j+1)));
                            elseif(j>floor(lwindowTimes/2))
                                tempSpikes=spikeTimes(spikeTimes>=windowTimes(j));
                                data(j)=sum(tempSpikes<windowTimes(j+1));
                        %                                     data(j) = sum((spikeTimes>=windowTimes(j) & spikeTimes<windowTimes(j+1)));
                            else
                                tempSpikes=spikeTimes(spikeTimes<windowTimes(j+1));
                                data(j)=sum(tempSpikes>=windowTimes(j));
                            end
                        end
                    
%                         timeVec    = round(timeVec*nstObj.sampleRate)/nstObj.sampleRate;
%                         for j=1:length(timeVec)-1 %number of bins
%                             if(j==(length(timeVec)-1))
%                                     data(j) = sum((spikeTimes>=timeVec(j) & spikeTimes<=timeVec(j+1)));
%                             else
%                                     data(j) = sum((spikeTimes>=timeVec(j) & spikeTimes<timeVec(j+1)));
%                             end
%                         end
%                          tV=repmat(timeVec,[length(spikeTimes) 1]);
%                          sT=repmat(spikeTimes',[1 length(timeVec)-1]);
%                          data= sT>=tV(:,1:end-1) & sT<tV(:,2:end);
%                          data(:,end)= sT(:,end)>=tV(:,end-1) & sT(:,end)<=tV(:,end);
%                          data=[sum(data) 0];
                        
                        sigRep = SignalObj(timeVec, data',nstObj.name,nstObj.xlabelval, nstObj.xunits, nstObj.yunits, nstObj.dataLabels);
                        nstObj.sigRep = sigRep;
                        if(max(sigRep.data)>1)
                            nstObj.isSigRepBin=0;
                        else
                            nstObj.isSigRepBin=1;
                        end
                    end
                    nstObj.sigRep = sigRep;
                    if(max(sigRep.data)>1)
                        nstObj.isSigRepBin=0;
                    else
                        nstObj.isSigRepBin=1;
                    end
            %end
        end
        function maxBinSize=getMaxBinSizeBinary(nstObj)
            % maxBinSize=getMaxBinSizeBinary(nstObj)
            % returns the maximum binsize or binwidth at which the
            % nspikeTrain still has a binary SignalObj representation
            if(length(nstObj.spikeTimes)>1)
                maxBinSize=min(diff(nstObj.spikeTimes));
            else
                maxBinSize=inf;
            end
        end
        function h = plotISISpectrumFunction(nstObj)
            figure;
            ISI=nstObj.getISIs;
            spikeTimes =nstObj.spikeTimes;
            h=plot(spikeTimes(2:end),ISI,'.');
            xlabel('time [s]');
            
            ylabel('ISI [s]');
            sigma=std(ISI);
            mu = mean(ISI);
            r= sigma/mu;
            B = (r-1)/(r+1);
            n=length(spikeTimes);
            An=(sqrt(n+2)*r-sqrt(n))./((sqrt(n+2)-2)*r+sqrt(n));
        end
        
        function windowedSpikeTimes = getSpikeTimes(nstObj, minTime,maxTime)
            if(nargin<3)
                maxTime = nstObj.maxTime;
            end
            if(nargin<2)
                minTime = nstObj.minTime;
            end
            index = and((nstObj.spikeTimes>=minTime),(nstObj.spikeTimes<=maxTime));
            windowedSpikeTimes = nstObj.spikeTimes(index);
        end
        function h = plotJointISIHistogram(nstObj)
            % based on: Detection of bursts in neuronal spike trains by the mean inter-spike interval method. (n.d.). Detection of bursts in neuronal spike trains by the mean inter-spike interval method.
            
            ISIs = nstObj.getISIs;
            meanISI = mean(ISIs);
            Ln = ISIs(ISIs<meanISI);
            ML = mean(Ln);
            loglog(ISIs(1:end-1),ISIs(2:end),'.'); hold on;
            v=axis;
            loglog(ML*[1;1],v(3:4),'k--');
            loglog(v(1:2),ML*[1;1],'k--');
            xlabel('ISI(t) [s]'); ylabel('ISI(t+1) [s]');
           
        end
        function fieldVal = getFieldVal(nstObj,fieldName)
            if(any(strcmp(fieldnames(nstObj),fieldName)))
                fieldVal = nstObj.(fieldName);
            else
                fieldVal = [];
            end
            
        end
        function counts = plotISIHistogram(nstObj,minTime,maxTime,numBins,handle)
%             if(nargin<6 || isempty(color))
%                 color = [0.831372559070587 0.815686285495758 0.7843137383461];
%             end
            if(nargin<5 || isempty(handle))
              handle=gca;
            end 
            if(nargin<4 || isempty(numBins))
                numBins = 100;
            end
            if(nargin<3 || isempty(maxTime))
                maxTime = nstObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime = nstObj.minTime;
            end
            
            
            ISIs = nstObj.getISIs;
%             index=and(ISIs>=minTime, ISIs<=maxTime);
%             ISIs = ISIs(index);
            binWidth=.001; %max(ISIs)/numBins;
%             binWidth=1/numBins;
            if(~isempty(ISIs))
                bins=0:binWidth:max(ISIs);

                %Make the ISI Histogram            
                counts = histc(ISIs,bins);

                %set(gcf,'CurrentAxes',handle);
                %bar(bins,histc(ISIs,bins)./sum(binWidth*counts),'histc');
                h=bar(bins,histc(ISIs,bins),'histc');
                set(h,'MarkerEdgeColor',[0 0 0],...
                'LineWidth',2,...
                'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461]);
            end
            
            
            axis tight;
            hx=xlabel('ISI [sec]');
%             hy=ylabel([nstObj.name ' counts']);
            hy=ylabel('Spike Counts');
            set([hx, hy],'FontName', 'Arial','FontSize',16,'FontWeight','bold');
            v=axis;
            axis tight;
%             axis([minTime, maxTime, v(3:4)]);
            
            %histfit(ISIs,numBins,'exponential');
%             h = get(gca,'Children');
%             set(h,'FaceColor',[.8 .8 1])

            
            %Fit exponential distribution to the data
%             [muhat,muci] = expfit(ISIs);
%             x=linspace(0,max(bins),1000);
%             y=exppdf(x,muhat);%;*1/muhat;
%             ci(:,1) = exppdf(x,muci(1));%*1/muci(1);
%             ci(:,2) = exppdf(x,muci(2));%*1/muci(2);
%             xpatch=[x fliplr(x)];
%             ypatch=[ci(:,1)',fliplr(ci(:,2)')];
            
%             hold on;
%             hfit1=plot(x,y,'r','Linewidth',3); %plot fit
%             p=patch(xpatch,ypatch,'r');
%             set(p,'facecolor','r','edgecolor','none'); %plot CI
%             alpha(.5);
            
            %Fit Weibul Distribution to Data
%             [parmhat,parmci] = wblfit(ISIs);
%             y = wblpdf(x,parmhat(1),parmhat(2));
%             hfit2=plot(x,y,'k','Linewidth',3); %plot fit
%             
            
            %Fit Gamma Distribution
%             [phat,pci] = gamfit(ISIs);
%             y = gampdf(x,phat(1),phat(2));
%             hfit3=plot(x,y,'g','Linewidth',3); %plot fit
%             
%             [phat,pci]=raylfit(ISIs)
%             y = raylpdf(x,phat) 
%         
%             numDigits=3;
%             expStr = ['exp: \lambda=' num2str(1/muhat,numDigits) ' [' num2str(1/muci(2),numDigits) ', ' num2str(1/muci(1),numDigits) ']'];
%              expStr = ['\lambda=' num2str(1/muhat,numDigits) ' [' num2str(1/muci(2),numDigits) ', ' num2str(1/muci(1),numDigits) ']'];
%             wblStr = ['weib: shape=' num2str(parmhat(2),numDigits) ', scale=' num2str(parmhat(1),numDigits)]; %scale and shape
%             gammaStr=['\Gamma: shape=' num2str(phat(1),numDigits) ', scale =' num2str(phat(2),numDigits)]; %shape and scale
%             legend([hfit1(1) hfit2(1) hfit3(1)] ,{expStr,wblStr,gammaStr});
%             
%             legend(hfit1(1) ,expStr);


%             subplot(1,2,2); nstObj.plotProbPlot(minTime,maxTime);
            
        end
        
        
        function h = plotExponentialFit(nstObj,minTime,maxTime,numBins, handle)
            if(nargin<5 || isempty(handle))
              handle=gca;
            end 
            if(nargin<4 || isempty(numBins))
                numBins = 200;
            end
            if(nargin<3 || isempty(maxTime))
                maxTime = nstObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime = nstObj.minTime;
            end
            h=figure;
            subplot(1,2,1); nstObj.plotISIHistogram(minTime,maxTime,numBins,handle)
            subplot(1,2,2); nstObj.plotProbPlot(minTime,maxTime,handle);
            
        end
        
        function h = plotProbPlot(nstObj,minTime,maxTime,handle)
            if(nargin<4 || isempty(handle))
                handle =gca;
            end
            if(nargin<3 || isempty(maxTime))
                maxTime = nstObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime = nstObj.minTime;
            end
            %set(gcf,'CurrentAxes',handle);
            ISIs = nstObj.getISIs(minTime,maxTime);
            h=probplot('exponential',ISIs);
            [muhat,muci] = expfit(ISIs);
%             hold on;
%             Z=1/muhat*(nstObj.spikeTimes(2:end) - nstObj.spikeTimes(1:end-1));
%             U = 1-exp(-Z); %store the rescaled spike times
% 
% 
%             KSSorted = sort( U,'ascend' );
%             N = length(KSSorted);
%             if(N~=0)
%                 xAxis=(([1:N]-.5)/N)';
%                 ks_stat = max(abs(KSSorted' - (([1:N]-.5)/N)')); 
%             else
%                 ks_stat=1;
%                 xAxis=[];
%             end
%             
%                          
% %             handle=plot(xAxis,KSSorted, 0:.01:1,0:.01:1, 'k-.',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
% 
%             %set(gca,'xtick',[],'ytick',[],'ztick', [])
% %             axis( [0 1 0 1] );
%             
%             xlabel('Uniform CDF');
%             ylabel('Empirical CDF of Rescaled ISIs');
%             title('KS Plot with 95% Confidence Intervals');
%             
            
            
            
        end
        function ISIs = getISIs(nstObj,minTime,maxTime)
            if(nargin<3 || isempty(maxTime))
                maxTime = nstObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime = nstObj.minTime;
            end
            spikeTimes = nstObj.getSpikeTimes(minTime,maxTime);
            ISIs = diff(spikeTimes)';
        
        end
        
        function minISI = getMinISI(nstObj,minTime, maxTime)
            if(nargin<3 || isempty(maxTime))
                maxTime = nstObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime = nstObj.minTime;
            end
            minISI = min(nstObj.getISIs(minTime,maxTime));
        end
        function nstCollObj = partitionNST(nstObj, windowTimes,normalizeTime,lbound,ubound)
            if(nargin<5 || isempty(ubound))
                if(nargin>4)
                    ubound = lbound;
                else
                    ubound=[];
                end
                
            end
            if(nargin<4 || isempty(lbound));
                lbound=[];
            end
            if(nargin<3 || isempty(normalizeTime))
                normalizeTime = [];
            end
            
%             nst = cell(1,length(windowTimes)-1);
            nst={};
            for i=1:length(windowTimes)-1
                
                minTime = round(windowTimes(i)*nstObj.sampleRate)/nstObj.sampleRate;
                maxTime = round(windowTimes(i+1)*nstObj.sampleRate)/nstObj.sampleRate;
%                 spikeTimes = round(nstObj.spikeTimes*nstObj.sampleRate)/nstObj.sampleRate;
                spikeTimes = nstObj.spikeTimes;%*nstObj.sampleRate)/nstObj.sampleRate;
                if(and(~isempty(lbound),~isempty(ubound)))
                    if(and(abs(maxTime-minTime)<=(ubound),abs(maxTime-minTime)>=(lbound)))
                        dim = length(nst);
                        if(i==(length(windowTimes)-1))
                             
                            spikeTimesSubset = spikeTimes(spikeTimes>=minTime & spikeTimes<=maxTime);
                        else
                            spikeTimesSubset = spikeTimes(spikeTimes>=minTime & spikeTimes<maxTime);
                        end
                        spikeTimesSubset = spikeTimesSubset - minTime;
%                         spikeTimes = round(spikeTimes*nstObj.sampleRate)/nstObj.sampleRate;
%                         minTime= round(minTime*nstObj.sampleRate)/nstObj.sampleRate;
%                         maxTime= round(maxTime*nstObj.sampleRate)/nstObj.sampleRate;
%                         timeVec    = round(timeVec*nstObj.sampleRate)/nstObj.sampleRate;
                        if(normalizeTime==1)
%                             nst{dim+1} = nspikeTrain(spikeTimes/max(spikeTimes),strcat([nstObj.name ',w' num2str(i)]));
                              spikeTimesSubset = spikeTimesSubset/(maxTime-minTime);
                              nst{dim+1} = nspikeTrain(spikeTimesSubset,nstObj.name);

                        else
%                             nst{dim+1} = nspikeTrain(spikeTimes,strcat([nstObj.name ',w' num2str(i)]));
                            nst{dim+1} = nspikeTrain(spikeTimesSubset,nstObj.name);
                        end
                    end
                    
                else
                     if(i==(length(windowTimes)-1))
                             
                        spikeTimesSubset = spikeTimes(and(spikeTimes>=minTime,spikeTimes<=maxTime));
                        spikeInInterval{i}= spikeTimesSubset;
                     else
                        spikeTimesSubset = spikeTimes(and(spikeTimes>=minTime,spikeTimes<maxTime));
                         spikeInInterval{i}= spikeTimesSubset;
                     end
                        
%                     spikeTimes = spikeTimes(and(spikeTimes>minTime,spikeTimes<=maxTime));
                    spikeTimesSubset = spikeTimesSubset - minTime;
%                     spikeTimes = round(spikeTimes*nstObj.sampleRate)/nstObj.sampleRate;
%                     %spikeTimes = spikeTimes./(windowTimes(i+1)-windowTimes(i));
%                     minTime= round(minTime*nstObj.sampleRate)/nstObj.sampleRate;
%                     maxTime= round(maxTime*nstObj.sampleRate)/nstObj.sampleRate;
                    if(normalizeTime==1)
%                         nst{i} = nspikeTrain(spikeTimes/max(spikeTimes),strcat([nstObj.name ',w' num2str(i)]));
%                         spikeTimes=round(spikeTimes/(maxTime-minTime)*nstObj.sampleRate)/nstObj.sampleRate;
                        spikeTimesSubset = spikeTimesSubset/(maxTime-minTime);
                        nst{i} = nspikeTrain(spikeTimesSubset,nstObj.name);
                    else
%                         nst{i} = nspikeTrain(spikeTimes,strcat([nstObj.name ',w' num2str(i)]));
                        nst{i} = nspikeTrain(spikeTimesSubset,nstObj.name);
                    end
                end
            end
            nstCollObj = nstColl(nst);
            if(normalizeTime==1)
                nstCollObj.setMinTime(0);
                nstCollObj.setMaxTime(1);
            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function answer = isSigRepBinary(nstObj)
            % answer = isSigRepBinary(nstObj) 
            % answer = 1 if the SignalObj representation is binary
            % answer = 0 otherwise
%             if(~isempty(nstObj.sigRep))
%               if(max(nstObj.sigRep.data)>1)
%               %if(max(nstObj.getSigRep.data)>1)
%                   answer=0;
%               else
%                   answer=1;
%               end      
%             else
              if(isempty(nstObj.sigRep) || (nstObj.sampleRate~=nstObj.sigRep.sampleRate))
                  nstObj.getSigRep;
              end
              if(max(nstObj.sigRep.data)>1)
                  answer=0;
              else
                  answer=1;
              end      
%             end
        end
        function rateSignal = computeRate(nstObj)
            % rateSignal = computeRate(nstObj)
            % generate rate function signal for the corresponding spikeTimes
            % not yet implemented
        end
        function restoreToOriginal(nstObj)
            % restoreToOriginal(nstObj)
            % restores the signalRep of the nspikeTrain to its original
            % state. Sets sampleRate to the original sampleRate of the
            % SignalObj. 
            % 
            % nstObj.minTime = min(nstObj.spikeTimes)
            % nstObj.maxTime = max(nstObj.spikeTimes)
            %
            %nstObj.sigRep.restoreToOriginal;
            %nstObj.isSigRepBin=nstObj.isSigRepBinary;
            %nstObj.sampleRate = nstObj.sigRep.sampleRate;
            nstObj.minTime = min(nstObj.spikeTimes);
            nstObj.maxTime = max(nstObj.spikeTimes);
        end
        function nCopy = nstCopy(nstObj)
           % nCopy = nstCopy(nstObj)
           % nCopy is a copy of the nspikeTrain nstObj. This function is
           % important since nspikeTrains have handle behavior. For
           % example,
           % n2= n1; %where n1 is a nspikeTrain
           % n2.setMinTime(-10);% also sets the minTime of n1 to -10
           % because both reference the same data.
           %
           % To avoid this:
           % n2=n1.nstCopy;
           % n2.setMinTime(-10); %only changes n2.
           name       = nstObj.name;
           sampleRate = nstObj.sampleRate;
           spikeTimes = nstObj.spikeTimes;
           minTime    = nstObj.minTime;
           maxTime    = nstObj.maxTime;
%            if(~isempty(nstObj.sigRep))
%             sig        = nstObj.sigRep.copySignal;
%            else
%             sig        = nstObj.getSigRep;  
%            end
           nCopy.sigRep = nstObj.getSigRep; 
           sig = nstObj.getSigRep;
           xlabelval  = sig.xlabelval;
           xunits     = sig.xunits;
           nCopy=nspikeTrain(spikeTimes,name,1/sampleRate,minTime,maxTime, xlabelval,xunits);
%            nCopy.sigRep = nstObj.getSigRep; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h=plot(nstObj, dHeight, yOffset, currentHandle) 
            % h=plot(nstObj, dHeight, yOffset, currentHandle) 
            % only plots the actual spikeTimes regardless of the signal
            % representation
            % dHeight: height of a spike on the plot
            % yOffset: offset to be used when plotting this nspikeTrain
            %          nstColl uses the offset for drawing an entire
            %          collection of spikes that were recorded
            %          simultaneously.
            if (nargin < 4)
                currentHandle = gca;
            end                   
            if (nargin < 3)
                yOffset = .5;
            end
            if (nargin < 2)
                dHeight = 1; 
            end
             
            if(~isempty(nstObj.spikeTimes))
                time   =   ones(2,1)*nstObj.spikeTimes;
            else
                time   = [];
            end
            
            spikes =  ([-dHeight/2; +dHeight/2]+yOffset)*ones(1,length(nstObj.spikeTimes));
            %h = figure(currentHandle);
            hold on;
            h=plot(currentHandle,time,spikes,'k');
            
            if(nargin==1) %Only plot labels if being plotted separately
                sig       = nstObj.getSigRep;
                xlabelval = sig.xlabelval;
                xunits    = sig.xunits;
                name      = sig.name;
                yunits    = sig.yunits;
                if(~strcmp(xunits,''))
                    xunitsStr=strcat(' [',xunits,']');
                else
                    xunitsStr='';
                end
                xlabel(strcat(xlabelval,xunitsStr));

                if(~strcmp(yunits,''))
                    yunitsStr=strcat(' [',yunits,']');
                else
                    yunitsStr='';
                end
                ylabel(strcat(name,yunitsStr));
                v=axis;
                if(nstObj.minTime~=nstObj.maxTime)
                    axis([nstObj.minTime,nstObj.maxTime,v(3),v(4)]);
                end
            end

            
        end
        
        function structure = toStructure(nstObj)
            fnames = fieldnames(nstObj);
            
            for i=1:length(fnames)
               currObj = nstObj.(fnames{i});
               if(isa(currObj,'double') || isa(currObj,'char'))
                   structure.(fnames{i}) = currObj;
               elseif(isa(currObj,'SignalObj'))
                   structure.(fnames{i}) = currObj.dataToStructure;
               end
                
                
            end
        end
        
    end
    methods (Static)
        function nstObj = fromStructure(structure)
            spikeTimes = structure.spikeTimes;
            name       = structure.name;
            binwidth   = 1/structure.sampleRate;
            minTime    = structure.minTime;
            maxTime    = structure.maxTime;
            xlabelval  = structure.xlabelval;
            xunits     = structure.xunits;
            yunits     = structure.yunits;
            dataLabels = structure.dataLabels;
            nstObj=nspikeTrain(spikeTimes,name,binwidth,minTime,maxTime,xlabelval, xunits, yunits,dataLabels);
        end
    end
    
end

