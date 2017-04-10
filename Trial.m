classdef Trial <handle
% TRIAL A trial has handle behavior to save memory. It can be reset to its orginal
% state by calling restoreToOriginal;
%
% Usage:
%       tObj=Trial(nspikeColl, covarColl, event, hist)
%           hist:       History Object. Optional argument;
%           event:      Events Object. Optional argument;
%           covarColl:  CovColl Object. Required.
%           nspikeColl: nstColl Obj
%
% <a href="matlab: methods('Trial')">methods</a>
% <a href="matlab:web('TrialExamples.html', '-helpbrowser')">Trial Examples</a> 
%
% see also <a href="matlab:help('CovColl')">CovColl</a>, 
% <a href="matlab:help('Covariate')">Covariate</a>, 
% <a href="matlab:help('SignalObj')">SignalObj</a>,
% <a href="matlab:help('nspikeTrain')">nspikeTrain</a>
%
% Reference page in Help browser
% <a href="matlab: doc('Trial')">doc Trial</a>

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

    properties
        % A trial consists of 
        nspikeColl; % a collection of neural spikes
        covarColl;  % a collection covariates
        ensCovHist; % the history structure used to create the ensCovColl
        ev;         % events
        history ;   % current history settings for this trial
        sampleRate; % sampleRate for all the covariates and nstColl
        minTime;    % minTime of all data or minTime of interest
        maxTime;    % maxTime of all data or maxTime of interest
        covMask;    % mask indicating visible covariates
        ensCovMask; % mask indicating which neurons to consider for neighbor effcts;
        neuronMask; % mask indicating visible neurons
        trainingWindow   % start and endtimes for training data
        validationWindow % start and endtime for validation data
       
    end
    properties (Hidden)
        ensCovColl; % a covariate collection of neighboring neuron spiking activity
    end
    

    methods
        function tObj=Trial(nspikeColl, covarColl, event, hist,ensCovHist,ensCovMask)
            % tObj=Trial(nspikeColl, covarColl, event, hist,ensCovHist)
            % nspikeColl: is an <a href="matlab:help('nstColl')">nstColl</a> object 
            %             containing all of the relevant spike trains for this experimental trial.
            % covarColl: is a <a href="matlab:help('CovColl')">CovColl</a>
            %            object containing all of the covariates associated with the trial
            % event:     <a href="matlab:help('Events')">Events</a> object. 
            % hist:      <a href="matlab:help('History')">History</a>
            %            object for any given spike train
            % ensCovHist:<a href="matlab:help('History')">History</a> 
            %            object that specifies how much history to include for the ensemble effect. 
            if(nargin<6)
                nSpikes = nspikeColl.numSpikeTrains;
                ensCovMask = ones(nSpikes,nSpikes)-eye(nSpikes,nSpikes);
            end
            if(nargin<5)
               ensCovHist = []; 
            end
            if(nargin<4)
                hist = [];
            end
            if(nargin<3)
                event=[];
            end
            
            if(isa(nspikeColl,'nstColl'))
                tObj.nspikeColl = nspikeColl;
            else
                error('nstColl is a required argument');
            end
            
            if(isa(covarColl,'CovColl'))
                tObj.covarColl  = covarColl;
            else
                error('CovColl is a required argument');
            end
            
            tObj.setTrialEvents(event);
            
    
            
            if(isa(hist,'History') || isa(hist,'double'))
                tObj.setHistory(hist);
            else
                tObj.history = [];
            end
            
            if(isa(ensCovHist,'History')|| isa(ensCovHist,'double'))
                tObj.setEnsCovHist(ensCovHist);
            else
                tObj.ensCovHist=[];
            end

            tObj.covMask    = covarColl.covMask;
            tObj.neuronMask = nspikeColl.neuronMask;
            tObj.ensCovMask = ensCovMask;
             if(~tObj.isSampleRateConsistent)
                tObj.makeConsistentSampleRate;
             else
                 tObj.sampleRate = tObj.covarColl.sampleRate;
             end
            
            tObj.makeConsistentTime;
            tObj.setTrialPartition([]); %default to all training data
            tObj.setTrialTimesFor('training');
%             tObj.setBatchMode('off'); % Turn batchMode off by default
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setTrialEvents(tObj,event)
            if(isa(event,'Events'))
                tObj.ev = event;
            else
                tObj.ev = [];
            end
        end
%         
%         function setBatchMode(tObj,mode)
%            if(isa(mode,'char'))
%                if(strcmp(mode,'on'))
%                 tObj.batchMode = 'on';
%                elseif(strcmp(mode,'off'));
%                 tObj.batchMode ='off';
%                end
%            elseif(isa(mode,'double'))
%                if(mode==1)
%                    tObj.batchMode='on';
%                elseif(mode==0)
%                    tObj.batchMode='off';
%                end
%            end
%                          
%         end
%         
        function setTrialPartition(tObj,partitionTimes)
            %partitionTimes is an array containing three or four elements
            % 4 elements :[startTraining stopTraining startValidation endValidation]
            % 3 elements [startTraining stopTraining/startValidation endValidation]
            % the middle value is used for end of training period and start
            % of validation period
            if(nargin<2 || isempty(partitionTimes))
                partitionTimes = tObj.getTrialPartition;
            end
            
            if(~isempty(partitionTimes))
                %partitionTimes = round(partitionTimes*tObj.sampleRate/10)./tObj.sampleRate*10; %make sure that the partition time are achievable at the current sampleRate;
                % if sampleRate is 100 then we multiply partionTimes by 10
                % and so we keep only the first decimal place.
                if(length(partitionTimes)==4)
                    trainingWindow   = [partitionTimes(1) partitionTimes(2)]; 
                    validationWindow = [partitionTimes(3) partitionTimes(4)]; 
                elseif(length(partitionTimes)==3)
                    trainingWindow   = [partitionTimes(1) partitionTimes(2)]; 
                    validationWindow = [partitionTimes(2) partitionTimes(3)]; 
                else
                    error('partitionTimes must be length 3 or 4');
                end

                tObj.trainingWindow   = trainingWindow;
                tObj.validationWindow = validationWindow;
                tObj.setMinTime(trainingWindow(1));
                tObj.setMaxTime(trainingWindow(2)); %default to ready for training
            end
        end
        function p=getTrialPartition(tObj)
            p1=tObj.trainingWindow;
            p2=tObj.validationWindow;
            p=[p1 p2];
            if(isempty(p))
                p1   = [tObj.minTime, tObj.maxTime];
                p2   = [tObj.maxTime, tObj.maxTime];
                p    = [p1, p2];
            end
        end
        
        
        function setTrialTimesFor(tObj,partitionName)
            if(nargin<2)
                partitionName = 'training';
            end
            
            p = tObj.getTrialPartition;
            
            if(strcmp(partitionName,'training'))
                timeWindow = p(1:2);
            elseif(strcmp(partitionName,'validation'))
                timeWindow = p(3:4);
            else
                error('partitionName must be either training or validation');
            end
            
            if(~isempty(timeWindow))
               %currSamplingRate = tObj.sampleRate;
               %tObj.restoreToOriginal;
               %tObj.resample(currSamplingRate); 
               tObj.setMinTime(timeWindow(1));
               tObj.setMaxTime(timeWindow(2));
               %  tObj.makeConsistentTime;
               %  tObj.restoreToOriginal;
               %  tObj.resample(currSamplingRate);                
            end
            
        end
        function setMinTime(tObj,minTime)
            % setMinTime(tObj,minTime) 
            % sets the minTime of interest of the trial to minTime
            
            if(nargin<2)
                minTime = tObj.findMinTime;
            end;
            tObj.nspikeColl.setMinTime(minTime);
            tObj.covarColl.setMinTime(minTime);
            if(~isempty(tObj.ensCovColl))
                tObj.ensCovColl.setMinTime(minTime);
            end
            %tObj.covarColl.covArray{1}.standardRep
            tObj.minTime = minTime;
         end
        function setMaxTime(tObj,maxTime)
            % setMaxTime(tObj,maxTime)
            % sets the maxTime of interest of the trial to maxTime
            if(nargin<2)
                maxTime=tObj.findMaxTime;
            end
            tObj.nspikeColl.setMaxTime(maxTime);
            tObj.covarColl.setMaxTime(maxTime);
            if(~isempty(tObj.ensCovColl))
                tObj.ensCovColl.setMaxTime(maxTime);
            end
            tObj.maxTime = maxTime;
        end
        
        function updateTimePartitions(tObj)
            if((~isempty(tObj.minTime))&& (~isempty(tObj.maxTime))) %avoid calling before maxTime and minTime are set
                p = tObj.getTrialPartition;
                if(~isempty(p))
                    training = p(1:2);
                    validation=p(3:end);
                    newTrainMin = max(tObj.minTime,training(1));
                    newTrainMax = min(tObj.maxTime,training(2));
                    newValMin   = max(tObj.minTime,validation(1));
                    newValMax   = min(tObj.maxTime,validation(2));

                    tObj.setTrialPartition([newTrainMin newTrainMax newValMin newValMax]);
                end
            end
            
        end
        
        function setSampleRate(tObj,sampleRate)
            % setSampleRate(tObj,sampleRate)
            % sets the sampleRate of the trial and all of its components to
            % sampleRate
              tObj.sampleRate=sampleRate;
              tObj.nspikeColl.resample(sampleRate);
              tObj.covarColl.resample(sampleRate);
              tObj.resampleEnsColl;
        end
        function setEnsCovMask(tObj,mask)
           % setEnsCovMask(tObj, mask)
           % sets the mask of neighboring neurons to be considered when
           % ensCovHist is set
           if(nargin<2 || isempty(mask))
               nSpikes = tObj.nspikeColl.numSpikeTrains;
               mask = ones(nSpikes,nSpikes)-eye(nSpikes,nSpikes);
           end
           tObj.ensCovMask=mask; % needs to be a nSpikeTrain x nSpikeTrain matrix with zeroes along diagonal if all neurons are possible neighbors
        end
        function setCovMask(tObj,mask)
            % setCovMask(tObj,mask)
            % sets the covariate mask of the trial and of the covColl to
            % mask
            if(isa(mask,'char'))
                if(strcmp(mask,'all'));
                    tObj.covarColl.resetMask;
                end
            else
                tObj.covarColl.setMask(mask);
            end
                tObj.covMask = tObj.covarColl.covMask;          
           
        end
        function setNeuronMask(tObj,mask)
            % setNeuronMask(tObj,mask)
            % sets the neuron mask of the trial and of the nstColl to mask
            tObj.nspikeColl.setMask(mask);
            tObj.neuronMask = tObj.nspikeColl.neuronMask;
        end
        function setNeighbors(tObj,varargin)
            tObj.nspikeColl.setNeighbors(varargin{:});
        end
        function setHistory(tObj,hist)
            % setHistory(tObj,hist)
            % sets the history object of the trial to hist.
            % hist can be of class History or a vector of doubles
            % specifying the windowTimes for the History object
            if(isempty(hist))
                tObj.history = [];
            else
                if(isa(hist,'History'))
                    tObj.history = hist;
                elseif(isa(hist,'cell'))
                    if(isa(hist{1},'History'))
                        tObj.history = hist; 
                    end
                elseif(isa(hist,'double')) %then we got windowTimes
                    [l,w]=size(hist);
                    if(min(l,w)>1)
                        error('Only one of the dimension of the windowTimes can be greater than 1.');
                    else
                        if(length(hist)>1)
                            tObj.history = History(hist);
                        else
                            error('At least two times points must be specified to determine a window');
                        end

                    end
                else
                    error('Can only set trial history to be an instantiation of the History object class or by using windowTimes');
                end
            end
        end      
        function setEnsCovHist(tObj, hist)
            % ensCovHist(tObj,hist)
            % sets the ensCovHist of the trial to hist.
            % hist can be of class History or a vector of doubles
            % specifying the windowTimes for the History object
            if(nargin<2)
                hist =[];
            end
               
            if(isempty(hist))
                tObj.ensCovHist = [];
                tObj.ensCovColl = [];
            else
                if(isa(hist,'History'))
                    tObj.ensCovHist = hist;
                    
                elseif(isa(hist,'double')) %then we got windowTimes
                    [l,w]=size(hist);
                    if(min(l,w)>1)
                        error('Only one of the dimension of the windowTimes can be greater than 1.');
                    else
                        if(length(hist)>1)
                            tObj.ensCovHist = History(hist);
                            
                        else
                            error('At least two times points must be specified to determine a window');
                        end

                    end
                else
                    error('Can only set trial ensCovHist to be an instantiation of the History object class or by using windowTimes');
                end
                % getEnsembleNeuronCovariates(nstCollObj,neuronNum,neighbor
                % Index,windowTimes)
                tObj.ensCovColl=tObj.getEnsembleNeuronCovariates(1,[],tObj.ensCovHist);
                
                % initialize ensCovColl to first neuron and its default
                % neighbors
            end
        end        
        function answer=isNeuronMaskSet(tObj)
            % answer=isNeuronMaskSet(tObj)
            % 1 if neuronMask of nstColl is set;
            % 0 otherwise
            answer=tObj.nspikeColl.isNeuronMaskSet;
        end
        function answer=isCovMaskSet(tObj)
            % answer=isCovMaskSet(tObj)
            % 1 if the covMask of the covColl is set
            % 0 otherwise
            answer = tObj.covarColl.isCovMaskSet;
        end
        function answer=isMaskSet(tObj)
            % answer=isMaskSet(tObj)
            % 1 if either neuronMask or covMask is set for the trial
            % 0 if neither is set
            answer=tObj.isNeuronMaskSet||tObj.isCovMaskSet;
        end
        function answer=isHistSet(tObj)
            % answer=isHistSet(tObj)
            % 1 if the history object of the trial is non-empty
            % 0 if it is empty
            if(~isempty(tObj.history))
               if(isa(tObj.history,'History'))
                   answer =1;
               elseif(isa(tObj.history,'cell'))
                   if(isa(tObj.history{1},'History'))
                       answer=1;
                   else
                       answer=0;
                   end
               end
            else
               answer =0; 
            end
            
            
        end
        function answer=isEnsCovHistSet(tObj)            
            % answer=isHistSet(tObj)
            % 1 if the history object of the trial is non-empty
            % 0 if it is empty
            answer=(~isempty(tObj.ensCovHist))&& (isa(tObj.ensCovHist,'History'));
        end
        
        function addCov(tObj,cov)
            tObj.covarColl.addToColl(cov);
            tObj.covMask = tObj.covarColl.covMask;
            if(~tObj.isSampleRateConsistent)
                tObj.makeConsistentSampleRate;
            end
            tObj.makeConsistentTime;
        end
        function removeCov(tObj,identifier)
            tObj.covarColl.removeCovariate(identifier);
            tObj.covMask = tObj.covarColl.covMask;
            if(~tObj.isSampleRateConsistent)
                tObj.makeConsistentSampleRate;
            end
            tObj.makeConsistentTime;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y=getSpikeVector(tObj,varargin)
            % y=getSpikeVector(tObj,varargin)
            % returns a matrix of spikes by calling the dataToMatrix method
            % of nstColl with any arguments passed in
            % see also <a href="matlab:help('nstColl.dataToMatrix')">nstColl.dataToMatrix</a>,
            y=tObj.nspikeColl.dataToMatrix(varargin{:});
        end   
        function X=getDesignMatrix(tObj,neuronNum,varargin)
            % X=getDesignMatrix(tObj,varargin)
            % returns a matrix X by calling covColl.getDataMatrix with repType='zero-mean'
            % along with any parameters passed in. 
            %
            % dataSelector must be in the following format
            %   dataSelector{1} = {'Position','x','y'};
            %   dataSelector{2} = {'Force','fx','fy','fz'};
            % see also <a
            % href="matlab:help('covColl.dataToMatrix')">covColl.dataToMatrix</a>
            
            if(nargin<2)
                error('Must specify neuronNumber to getDesignMatrix!');
            end
                    
            %repType='zero-mean';
            repType='standard';
            X=tObj.covarColl.dataToMatrix(repType,varargin{:});
            
            if(tObj.isHistSet)
                H=tObj.getHistMatrices(neuronNum);
                if(isempty(X))
                    X=H;
                else
                    X=[X,H];
                end
            end
            if(tObj.isEnsCovHistSet)
                includedNeurons = find(tObj.ensCovMask(:,neuronNum)==1);
                E = tObj.getEnsCovMatrix(neuronNum,includedNeurons);
                if(isempty(X))
                    X=E;
                else
                    X=[X,E];
                end
            end
        end        
        function ensCovMatOut = getEnsCovMatrix(tObj,neuronNum,includedNeurons,varargin)
            if(nargin<3)
               includedNeurons = find(tObj.ensCovMask(:,neuronNum)==1);
            end
            if(nargin<2)
                error('Must specify neuronNumber to get the right neighbors!');
            end
            if(tObj.isEnsCovHistSet && ~isempty(neuronNum))
                if(neuronNum>0 && neuronNum<=tObj.nspikeColl.numSpikeTrains)
                    ensCovCollTemp = tObj.ensCovColl;
                    neighbors = tObj.getNeuronNeighbors(neuronNum);
                    neighbors = intersect(neighbors, includedNeurons);
                    ensCovCollTemp.maskAwayAllExcept(neighbors);
                    repType='standard';
                    ensCovMatOut=ensCovCollTemp.dataToMatrix(repType,varargin{:});
                end
            else
                ensCovMatOut=[]; %dont return matrix if neuron number not specified.
               % display('isEnsCovEnabled=0 - empty matrix returnecd');
            end
        end
        function histCovColl = getHistForNeurons(tObj,neuronIndex)
            % histCovColl = getHistForNeurons(tObj,neuronIndex)
            % returns a CovColl with number of Covariates equal to the
            % length of neuronIndex.
            % CovColl.getCov(i) corresponds to the History for the
            % nspikeTrain corresponding to neuronIndex(i)
            %
            if(tObj.isHistSet) 
                nst=tObj.nspikeColl.getNST(neuronIndex);
                if(length(neuronIndex)>=1)
                    if(isa(tObj.history,'History'))
                        histCovColl = tObj.history.computeHistory(nst);
                    elseif(isa(tObj.history,'cell'))
                        for i=1:length(tObj.history)
                            if(i==1)
                                histCovColl = tObj.history{1}.computeHistory(nst,i);
%                                 histCovColl.getCov
                            else
                                tempHistCovColl = tObj.history{i}.computeHistory(nst,i);
                                histCovColl.addToColl(tempHistCovColl);
                            end
                        end
                    end
                    
                end
            else
                histCovColl = []; %returns an empty array to indicate no History
                display('History has not been specified. Empty array returned');
                display('Set Trial history and retry');
            end
            if(~isempty(histCovColl))
                if(isa(nst,'cell'))
                    histCovColl.setSampleRate(nst{1}.sampleRate);
                else
                    histCovColl.setSampleRate(nst.sampleRate);
                end
            end
            
        end
        function matrices = getHistMatrices(tObj,neuronIndex)
           % matrices = getHistMatrices(tObj,neuronIndex)
           % return a matrix representation of the History Covariates
           % corresponding the History object specified for the trial and
           % the data corresponding to neuronIndex
           if(tObj.isHistSet)
               histCovColl = tObj.getHistForNeurons(neuronIndex);
               matrices = cell(1,length(neuronIndex));
               if(length(histCovColl)==1)
                   matrices = histCovColl.dataToMatrix;
               else
                   for i=1:length(histCovColl)
                        matrices{i} = histCovColl{i}.dataToMatrix;
                   end
               end
           else
               matrices = cell(1, length(neuronIndex));
               for i=1:length(neuronIndex)
                    matrices{i} = zeros(length(tObj.nspikeColl.getNST(neuronIndex(i)).sigRep.time),0);
               end
           end
        end
        function ensembleCovariates = getEnsembleNeuronCovariates(tObj,varargin)
            % getEnsembleNeuronCovariates(nstCollObj,neuronNum,neighborIndex,windowTimes)
                 ensembleCovariates = tObj.nspikeColl.getEnsembleNeuronCovariates(varargin{:});
        end
        function index    = getNeuronIndFromMask(tObj)
            % index    = getNeuronIndFromMask(tObj)
            % see also <a href="matlab:help('nstColl.getIndFromMask')">nstColl.getIndFromMask</a>
            index=tObj.nspikeColl.getIndFromMask;
        end
        function num = getNumUniqueNeurons(tObj)
           num = length(tObj.nspikeColl.uniqueNeuronNames);
        end
        function names = getNeuronNames(tObj)
            names = tObj.nspikeColl.getNSTnames;
        end
        function unames = getUniqueNeuronNames(tObj)
            unames = tObj.nspikeColl.getUniqueNSTnames;
        end
        function index = getNeuronIndFromName(tObj,neuronName)
            tempInd = tObj.nspikeColl.getNSTIndicesFromName(neuronName);
            currMask = find(tObj.neuronMask==1);
            if(isa(tempInd,'cell'))
                for i=1:length(tempInd)
                    index{i} = intersect(currMask,tempInd{i});
                end
            elseif(isa(tempInd,'double'))
                index = intersect(currMask,tempInd);
            end
            
                
        end
        
       
        function n        = getNeuronNeighbors(tObj, neuronNum)
            if(nargin<2)
                neuronNum = tObj.getNeuronIndFromMask;
            end
            n=tObj.nspikeColl.getNeighbors(neuronNum);
        end
        
        function selector = getCovSelectorFromMask(tObj)
            % selector = getCovSelectorFromMask(tObj)
            % see also <a href="matlab:help('CovColl.getSelectorFromMasks')">CovColl.getSelectorFromMasks</a>
            selector = tObj.covarColl.getSelectorFromMasks;
        end
        function cov      = getCov(tObj,identifier)
            % cov = getCov(tObj,identifier)
            % see also <a href="matlab:help('CovColl.getCov')">CovColl.getCov</a>
            cov=tObj.covarColl.getCov(identifier);
        end
        function NST      = getNeuron(tObj,identifier)
            % NST = getNeuron(tObj,identifier)
            % see also <a href="matlab:help('nstColl.getNST')">nstColl.getNST</a>
            NST = tObj.nspikeColl.getNST(identifier);
        end
        function e        = getEvents(tObj)
            % e = getEvents(tObj)
            % Returns the Events object associated with the Trial
            % e is either [] or an Events object
            e=tObj.ev;
        end
        
        function l = getAllLabels(tObj)
            l=tObj.getAllCovLabels;
            offset=length(l);
            if(tObj.isHistSet)
                l2=tObj.getHistLabels;
                for i=1:length(l2)
                    l{offset+i} = l2{i};
                end
            end
            offset=length(l);
            if(tObj.isEnsCovHistSet)
                l3=tObj.getEnsCovLabels;
                for i=1:length(l3)
                    l{offset+i} = l3{i};
                end
            end
        end            
        function n = getNumHist(tObj)
            % n = getNumHist(tObj)
            % Return the number of history windows (length windowTimes - 1)
           if(tObj.isHistSet)
               if(isa(tObj.history,'History'))
                n=length(tObj.history.windowTimes)-1;
               elseif(and(isa(tObj.history,'cell'),isa(tObj.history{1},'History')))
                   for i=1:length(tObj.history)
                       n(i) = length(tObj.history{i}.windowTimes)-1;
                   end
               end
           else
               n=0;
           end
           
        end
        function l = getAllCovLabels(tObj)
            l =tObj.covarColl.getAllCovLabels;
        end
        function l = getCovLabelsFromMask(tObj)
            l = tObj.covarColl.getCovLabelsFromMask;
        end
        
        function l = getHistLabels(tObj)
            if(tObj.isHistSet)
                histCovColl = tObj.getHistForNeurons(1);
                if(isa(histCovColl,'Covariate'))
                   l= histCovColl.dataLabels;
                 elseif(isa(histCovColl,'CovColl'))
                     l = histCovColl.getAllCovLabels;
                 else
                     error('histCovColl must be either a Covariate or a CovColl');
                end
            else
                l={};
            end
            
        end        
        function l = getEnsCovLabels(tObj)
           if(tObj.isEnsCovHistSet)
                ensCovCollTemp = tObj.ensCovColl;
                if(isa(ensCovCollTemp,'CovColl'))
%                     if(nargin<2)
%                         %l = ensCovCollTemp.getCovLabelsFromMask; %return labels from current mask
%                         getAllCovLabels
%                     else %set mask and return those labels
%                          ensCovCollTemp.maskAwayAllExcept(tObj.getNeuronNeighbors(neuronNum));
%                          l = ensCovCollTemp.getCovLabelsFromMask;
%                     end
                    l=ensCovCollTemp.getAllCovLabels;
                 else
                     error('ensCovColl must be either a CovColl');
                end
           else
               l={};
           end
            
        end
        function l = getEnsCovLabelsFromMask(tObj,neuronNum)
           if(nargin<2)
               error('Must specify neuron number!');
           end
           
           if(isa(neuronNum,'char')) 
              tempIndices = tObj.getNeuronIndFromName(neuronNum);
              if(length(tempIndices)>1)
                 neuronNum=tempIndices(1);
              else
                  neuronNum=tempIndices;
              end
           end
           if(tObj.isEnsCovHistSet)
                ensCovCollTemp = tObj.ensCovColl;
                if(isa(ensCovCollTemp,'CovColl'))
%                     if(nargin<2)
%                         %l = ensCovCollTemp.getCovLabelsFromMask; %return labels from current mask
%                         getAllCovLabels
%                     else %set mask and return those labels
%                          ensCovCollTemp.maskAwayAllExcept(tObj.getNeuronNeighbors(neuronNum));
%                          l = ensCovCollTemp.getCovLabelsFromMask;
%                     end
                     if(neuronNum>0 && neuronNum<=tObj.nspikeColl.numSpikeTrains)
                        ensCovCollTemp.maskAwayAllExcept(tObj.getNeuronNeighbors(neuronNum));
                        l=ensCovCollTemp.getCovLabelsFromMask;
                     else
                        error('NeuronNum is out of bounds!');
                     end
                 else
                     error('ensCovColl must be either a CovColl');
                end
           else
               l={};
           end
            
        end
        function l = getLabelsFromMask(tObj,neuronNum)
            if(nargin<2)
                error('To get right labels need to specify neuronNum'); %because of ensemble
            end
            l = tObj.getCovLabelsFromMask;
            offset=length(l);
            if(tObj.isHistSet)
                l2=tObj.getHistLabels;
                l((1:length(l2))+offset) = l2;
            end
            offset=length(l);
            if(tObj.isEnsCovHistSet)
                l3=tObj.getEnsCovLabelsFromMask(neuronNum);
                 l((1:length(l3))+offset) = l3;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function shiftTrial(tObj,deltaT)
%             tObj.covarColl.shiftCovariates(deltaT);
%             tObj.nspikeColl.shift(deltaT);
%             tObj.makeConsistentTime;
%             tObj.setTrialPartition([]); %default to all training data
%             %tObj.setTrialTimesFor('training');
%         end
%         
        function flatMask=flattenCovMask(tObj)
             flatMask = tObj.covarColl.flattenCovMask;
%              if(tObj.isHistSet)
%                  histCovColl = tObj.getHistForNeurons(1);
%                  if(isa(histCovColl,'Covariate'))
%                     flatMask = [flatMask histCovColl.dataMask];
%                  elseif(isa(histCovColl,'CovColl'))
%                      flatMask = [flatMask histCovColl.flattenCovMask];
%                  else
%                      error('histCovColl must be either a Covariate or a CovColl');
%                  end
%              end
             
        end        
        function flatMask=flattenMask(tObj)
             flatMask = tObj.flattenCovMask;
             if(tObj.isHistSet)
                 histCovColl = tObj.getHistForNeurons(1);
                 if(isa(histCovColl,'Covariate'))
                    flatMask = [flatMask histCovColl.dataMask];
                 elseif(isa(histCovColl,'CovColl'))
                     flatMask = [flatMask histCovColl.flattenCovMask];
                 else
                     error('histCovColl must be either a Covariate or a CovColl');
                 end
             end
             if(tObj.isEnsCovHistSet)
                 ensCovCollTemp = tObj.ensCovColl;
                 if(isa(ensCovCollTemp,'CovColl'))
                     flatMask = [flatMask ensCovCollTemp.flattenCovMask];
                 else
                     error('ensCovCollTemp must be either a CovColl');
                 end
             end
        end        
        function shiftCovariates(tObj, varargin)
            tObj.covarColl.setCovShift(varargin{:});
            tObj.makeConsistentTime;
        end
        
        function resetEnsCovMask(tObj)
            nSpikes = tObj.nspikeColl.numSpikeTrains;
            tObj.ensCovMask = ones(nSpikes,nSpikes)-eye(nSpikes,nSpikes);
        end
        function resetCovMask(tObj)
            % resetCovMask(tObj)
            % see also <a href="matlab:help('CovColl.resetMask')">CovColl.resetMask</a>
            tObj.covarColl.resetMask;
        end
        function resetNeuronMask(tObj)
            % resetNeuronMask(tObj)
            % see also <a href="matlab:help('nstColl.resetMask')">nstColl.resetMask</a>
            tObj.nspikeColl.resetMask;
        end
        function resetHistory(tObj)
            % resetHistory(tObj)
            % Sets the History object associated with this Trial equal to
            % the empty array [].
            tObj.history = [];
        end
%         function shiftCovs(tObj, covLag)
%             % shiftCovs(tObj, covLag)
%             % if covLag is not specified or isempty ([]) all covariates are
%             % href="matlab:help('CovColl.restoreToOriginal')">CovColl.restoreToOriginal</a>
%             % If specified, all covariates are shifted by covLag
%             % href="matlab:help('Covariate.shift')">Covariate.shift</a>
%             %
%             % The trial minTime and maxTime properties are updated to
%             % reflect any changes.
%             
%             if(nargin<2)
%                 covLag = [];
%             end
%             if(isempty(covLag))
%                 for i =1:tObj.covarColl.numCov
%                     tCov = tObj.covarColl.getCov(i);
%                     tCov.restoreToOriginal;
%                 end
%             else
%                 for i =1:tObj.covarColl.numCov
%                     tCov = tObj.covarColl.getCov(i);
%                     tCov.shift(covLag);
%                 end
%             end
%             tObj.makeConsistentTime;
%         end
        function resample(tObj,sampleRate)
            % resample(tObj,sampleRate)
            % calls <a
            % href="matlab:help('Trial.setSampleRate')">Trial.setSampleRate</a>
            tObj.setSampleRate(sampleRate);
        end        
        function resampleEnsColl(tObj)
            if(~isempty(tObj.ensCovColl) && ~isempty(tObj.ensCovHist))
                 tObj.ensCovColl=tObj.getEnsembleNeuronCovariates(1,[],tObj.ensCovHist);
            else
                tObj.setEnsCovHist; %set to empty;
            end
        end
        function restoreToOriginal(tObj)
            % restoreToOriginal(tObj)
            % calls <a
            % href="matlab:help('nstColl.restoreToOriginal')">nstColl.resto
            % reToOriginal</a>, <a
            % href="matlab:help('CovColl.restoreToOriginal')">CovColl.restoreToOriginal</a>
            % resets minTime and maxTime of the Trial according to the
            % changes that happened when restoring the nstColl and CovColl
            % to their original states.
              tObj.nspikeColl.restoreToOriginal;
              tObj.covarColl.restoreToOriginal;
              if(~tObj.isSampleRateConsistent)
                tObj.makeConsistentSampleRate;
              end
              tObj.resampleEnsColl; % compute at the new sampling rate
              tObj.makeConsistentTime;
        end
        function makeConsistentSampleRate(tObj)
            tObj.resample(tObj.findMaxSampleRate);
        end            
        function makeConsistentTime(tObj)
              tObj.setMinTime;
              tObj.setMaxTime;
              %tObj.setTrialPartition;
              %tObj.updateTimePartitions;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotRaster(tObj,handle)
            if(nargin<2)
                handle =gcf;
            end
        end
        function plotCovariates(tObj,handle)
            if(nargin<2)
                handle =gcf;
            end
            numCovars = tObj.covarColl.nActCovar; %accounts for masks
            figure(handle);
            if(numCovars==1)
                tObj.covarColl.plot; hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot;
                end
            elseif(numCovars==2)
                a1=subplot(1,2,1); 
                a2=subplot(1,2,2);
                tObj.covarColl.plot([a1,a2]); hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot([a1,a2]);
                end
            elseif(numCovars==3)
                a1=subplot(3,2,[1 3 5]);
                tObj.nspikeColl.plot;hold on;
                a2=subplot(3,2,2); a3=subplot(3,2,4); a4=subplot(3,2,6);
                tObj.covarColl.plot([a2,a3 a4]); hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot([a1,a2,a3,a4]);
                end
            else
                figure(handle);
                tObj.nspikeColl.plot;  hold on; 
                if(~isempty(tObj.ev))
                    tObj.ev.plot;
                end
                figure;
                tObj.covarColl.plot;  hold on; 
                if(~isempty(tObj.ev))
                    tObj.ev.plot;
                end
            end
        end
        
        function plot(tObj,handle)
            % plot(tObj,handle)
            % plots the Trial on the figure handle specified.
            % if handle is not specified, then handle = gcf;
            if(nargin<2)
                handle =gcf;
            end
            numCovars = tObj.covarColl.nActCovar; %accounts for masks
            figure(handle);
            if(numCovars==1)
                a1=subplot(2,2,[1 3]);
                tObj.nspikeColl.plot;hold on;
                a2=subplot(2,2,[2 4]); 
                tObj.covarColl.plot(a2); hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot([a1,a2]);
                end
            elseif(numCovars==2)
                a1=subplot(2,2,[1 3]);
                tObj.nspikeColl.plot;hold on;
                a2=subplot(2,2,2); a3=subplot(2,2,4);
                tObj.covarColl.plot([a2,a3]); hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot([a1,a2,a3]);
                end
            elseif(numCovars==3)
                a1=subplot(3,2,[1 3 5]);
                tObj.nspikeColl.plot;hold on;
                a2=subplot(3,2,2); a3=subplot(3,2,4); a4=subplot(3,2,6);
                tObj.covarColl.plot([a2,a3 a4]); hold on;
                if(~isempty(tObj.ev))
                    tObj.ev.plot([a1,a2,a3,a4]);
                end
            else
                figure(handle);
                tObj.nspikeColl.plot;  hold on; 
                if(~isempty(tObj.ev))
                    tObj.ev.plot;
                end
                figure;
                tObj.covarColl.plot;  hold on; 
                if(~isempty(tObj.ev))
                    tObj.ev.plot;
                end
            end
            %make sure events plot on every covariate plot
            %add events and covariate labels to the legend or mark events
            %with text on the screen.
        end
        function structure = toStructure(tObj)
            fNames =fieldnames(tObj);
            for i=1:length(fNames)
               currObj = tObj.(fNames{i});
               if(isa(currObj,'double')||isa(currObj,'cell'))
                   structure.(fNames{i}) = currObj;
               elseif(isa(currObj,'nstColl')||isa(currObj,'CovColl')||isa(currObj,'Events')||isa(currObj,'History'))
                   if(~isempty(currObj))
                    structure.(fNames{i}) = currObj.toStructure;
                   else
                    structure.(fNames{i}) = currObj;
                   end
               end
            end
        end
           
    end
    methods (Static)
        function tObj = fromStructure(structure)
            nspikeColl = nstColl.fromStructure(structure.nspikeColl);
            covarColl  = CovColl.fromStructure(structure.covarColl);
            ev         = Events.fromStructure(structure.ev);
            h          = History.fromStructure(structure.history);
            ensHist    = History.fromStructure(structure.ensCovHist);
            tObj = Trial(nspikeColl,covarColl,ev, h, ensHist);
           
            minTime    = structure.minTime;
            maxTime    = structure.maxTime;
            
            tObj.setMinTime(minTime);
            tObj.setMaxTime(maxTime);
            trainingW  = structure.trainingWindow;
            validationW= structure.validationWindow;
            
            tObj.setTrialPartition([trainingW validationW]);
        end
    end
    methods (Access = private)

        function answer = isSampleRateConsistent(tObj)
            answer = (tObj.findMinSampleRate == tObj.findMaxSampleRate);
        end
        function minTime = findMinTime(tObj)
            minTime=min([tObj.nspikeColl.minTime,tObj.covarColl.minTime]);
        end        
        function maxTime = findMaxTime(tObj)
            maxTime =max([tObj.nspikeColl.maxTime, tObj.covarColl.maxTime]);
        end        
        function minSampleRate = findMinSampleRate(tObj)
            minSampleRate = min([tObj.sampleRate, tObj.nspikeColl.sampleRate, tObj.covarColl.sampleRate]);
        end        
        function maxSampleRate = findMaxSampleRate(tObj)
            maxSampleRate = max([tObj.sampleRate, tObj.nspikeColl.sampleRate, tObj.covarColl.sampleRate]);
        end        

    end
    
end


   
