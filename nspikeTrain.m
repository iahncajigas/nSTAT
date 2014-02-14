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
    end
    
    methods
        function nst=nspikeTrain(spikeTimes,name,binwidth,minTime,maxTime, varargin)
            %varargin: xlabelval, xunits, yunits,dataLabels to Signal
            %constructor
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
            nst.setSigRep(binwidth, minTime, maxTime,varargin{:});

%             nst.sigRep = [];
%             nst.isSigRepBin=[];
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
        function setName(nstObj,name)
           % setName(nstObj,name) 
           % set the name after construction
           nstObj.name = name;
%            if(~isempty(nstObj.sigRep))
%                 nstObj.sigRep.setName(name);
%            end
        end
        function sigRep = setSigRep(nstObj, varargin)
            %sigRep = setSigRep(nstObj, varargin)
            %varargin: binwidth,minTime,maxTime,xlabelval, xunits,yunits,dataLabels
              nstObj.sigRep = nstObj.getSigRep(varargin{:});
              nstObj.isSigRepBin = nstObj.isSigRepBinary;
              nstObj.sampleRate = nstObj.sigRep.sampleRate;
              sigRep=nstObj.sigRep;
              nstObj.minTime = nstObj.sigRep.minTime;
              nstObj.maxTime = nstObj.sigRep.maxTime;
        end        
        function setMinTime(nstObj,minTime)
           % setMinTime(sObj,nstObj)
           % sets the minimun value of the time vector for the SignalObj representation of the nspikeTrain.
           %nstObj.sigRep.setMinTime(minTime);
           nstObj.minTime=minTime;
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
            %s=nstObj.sigRep;
            %binwidth=1/sampleRate;
            %nstObj.setSigRep(binwidth,nstObj.minTime,nstObj.maxTime,s.xlabelval,s.xunits);
            nstObj.sampleRate = sampleRate;
        end            
        function sigRep=getSigRep(nstObj,binwidth,minTime,maxTime,varargin)
            % sigRep=getSigRep(nstObj,binwidth,minTime,maxTime,varargin)
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
                end
                  if(and(~isempty(maxTime),~isempty(minTime)))
                    timeVec=linspace(minTime,maxTime,ceil((1/binwidth)*abs(maxTime-minTime)/binwidth)*binwidth+1); %scaling by binwidth to avoid roundoff error
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
                            spikeTimes = round(spikeTimes*nstObj.sampleRate*2)/(nstObj.sampleRate*2);
                            windowTimes    = round(windowTimes*nstObj.sampleRate*2)/(2*nstObj.sampleRate);
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
                         
                            sigRep = SignalObj(timeVec, data',nstObj.name,varargin{:});
                            nstObj.sigRep = sigRep;
                            nstObj.isSigRepBin=nstObj.isSigRepBinary;
                        end 
                   else
                        %rounding avoids comparison errors due to
                        %differences in non-significant digits
                        spikeTimes = nstObj.spikeTimes;
                        spikeTimes = round(spikeTimes*nstObj.sampleRate*2)/(nstObj.sampleRate*2);
                        windowTimes    = round(windowTimes*nstObj.sampleRate*2)/(2*nstObj.sampleRate);
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
                        
                        sigRep = SignalObj(timeVec, data',nstObj.name,varargin{:});
                        nstObj.sigRep = sigRep;
                        nstObj.isSigRepBin=nstObj.isSigRepBinary;
                   end
                  nstObj.sigRep = sigRep;
                  nstObj.isSigRepBin=nstObj.isSigRepBinary;
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
        function counts = plotISIHistogram(nstObj,minTime,maxTime,numBins,handle)
%             if(nargin<6 || isempty(color))
%                 color = [0.831372559070587 0.815686285495758 0.7843137383461];
%             end
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
            
            
            ISIs = nstObj.getISIs(minTime,maxTime);
            
            binWidth=max(ISIs)/numBins;
%             binWidth=1/numBins;
            bins=0:binWidth:max(ISIs);

            %Make the ISI Histogram            
            counts = histc(ISIs,bins);

            %set(gcf,'CurrentAxes',handle);
            bar(bins,histc(ISIs,bins)./sum(binWidth*counts),'histc');
            set(get(gca,'Children'),'MarkerEdgeColor',[0 0 0],...
                'LineWidth',2,...
                'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461]);
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
            axis tight;
            hx=xlabel('ISI [sec]');
%             hy=ylabel([nstObj.name ' counts']);
            hy=ylabel('Spike Counts');
            set([hx, hy],'FontName', 'Arial','FontSize',16,'FontWeight','bold');

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
            ISIs = diff(spikeTimes);
        
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
            nstObj=nspikeTrain(spikeTimes,name,binwidth,minTime,maxTime);
        end
    end
    
end

