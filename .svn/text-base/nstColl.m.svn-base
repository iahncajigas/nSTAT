classdef nstColl < handle
% NSTCOLL A collection of nspikeTrains
%   
% <a href="matlab: methods('nstColl')">methods</a>    
% <a href="matlab:web('nstCollExamples.html', '-helpbrowser')">nstColl Examples</a> 
%
% see also <a href="matlab:help('CovColl')">CovColl</a>, <a
% href="matlab:help('Covariate')">Covariate</a>, <a
% href="matlab:help('SignalObj')">SignalObj</a>,<a
% href="matlab:help('nspikeTrain')">nspikeTrain</a>
%
% Reference page in Help browser
% <a href="matlab: doc('nstColl')">doc nstColl</a>

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
        nstrain; %An array of nspikeTrain objects
        numSpikeTrains; % a running count of how many nspikeTrains are in object
        minTime %Time first spike occurs in the collection
        maxTime %TIme last spike occurs in the collection
        sampleRate
        neuronMask
        neuronNames
%         isSigRepBinary
        neighbors %the ith row specifies neighbors of ith neuron
    end
    properties (Dependent = true)
       
       uniqueNeuronNames
    end
    
    methods
        function nstCollObj=nstColl(nst)
        % nstCollObj=nstColl(nst)
        % nst is a cell array of nspikeTrains, a single nspikeTrains
        % or not specified. If not specified, an empty nstColl object is
        % created.
                nstCollObj.numSpikeTrains = 0;
                nstCollObj.minTime=inf;
                nstCollObj.maxTime=-inf;
                nstCollObj.neuronMask=[];
                nstCollObj.sampleRate=-inf;
                nstCollObj.neighbors = [];
%                 nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
            if(nargin<1) %Then we were called without any data
                nstCollObj.nstrain=cell(1);
            else
                nstCollObj.addToColl(nst);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function minTime = getFirstSpikeTime(nstCollObj)
            % minTime = getFirstSpikeTime(nstCollObj)
            % minTime is the time of the first spike from all of the
            % nspikeTrains in the collection
            minTime=nstCollObj.minTime;
        end        
        function maxTime = getLastSpikeTime(nstCollObj)
            % maxTime = getLstSpikeTime(nstCollObj)
            % maxTime is the time of the last spike from all of the
            % nspikeTrains in the collection
            maxTime=nstCollObj.maxTime;
        end
%         function nstCollObj = shift(nstCollObj,deltaT)
%             minTime=inf;
%             maxTime=-inf;
%            for i=1:nstCollObj.numSpikeTrains
%                 nstCollObj.getNST(i).shift(deltaT);
%                 minTime=min(nstCollObj.getNST(i).minTime,minTime);
%                 maxTime=max(nstCollObj.getNST(i).maxTime,maxTime);
%             end
%             nstCollObj.setMinTime(minTime);
%             nstCollObj.setMaxTime(maxTime);
%         end
        function answer  = getMaxBinSizeBinary(nstCollObj)
            % answer  = getMaxBinSizeBinary(nstCollObj)
            % answer is the binsize above which at least one nspikeTrain
            % will no longer have a binary representation
            %
            % Only unmasked ('visible') nspikeTrains are considered.
            if(nstCollObj.isNeuronMaskSet)
                selectorArray = nstCollObj.getIndFromMask;
            else
                selectorArray = 1:nstCollObj.numSpikeTrains;
            end
            val = zeros(1,length(selectorArray));
            for i=1:length(selectorArray)
%                 i
                val(i) = nstCollObj.getNST(selectorArray(i)).getMaxBinSizeBinary;
            end
            answer = min(val);
        end        
%         function value = get.isSigRepBinary(nstCollObj)
%             % value is 1 if all of the nspikeTrains in the collection have
%             % a binary representation. value=0 otherwise.
% %             value =nstCollObj.BinarySigRep;
%         end   
        function uniqueNames = get.uniqueNeuronNames(nstCollObj)
            uniqueNames=nstCollObj.getUniqueNSTnames;
        end

       
        function [n numNeighbors] = getNeighbors(nstCollObj, neuronNum)
           if(length(neuronNum)==1)
                %neuronNum: if neuronNum is not a scalar, a cell array is
                %returned. The ith entry of the cell array is a row vector
                %with the indicies of the ith neuron. If all of the rows
                %have the same length, then a matrix is returned instead of
                %a cell array.
               if(~nstCollObj.areNeighborsSet)
                    nstCollObj.setNeighbors; %default behavior
               end
                    availNeurons =nstCollObj.getIndFromMaskMinusOne(neuronNum);
                    if(isa(nstCollObj.neighbors,'cell'))
                        nTemp = nstCollObj.neighbors{neuronNum};
                    else
                        nTemp=nstCollObj.neighbors(neuronNum,:);
                    end
                    offset=0;
                    n=[];
                    for i=1:length(nTemp)
                        if(any(nTemp(i)==availNeurons))
                           offset=offset+1;
                           n(offset)=nTemp(i);
                        end
                    end
                    numNeighbors = length(n);
           else %if more than 1 neuronNum
               for i=1:length(neuronNum) % call above for each
                   [nTemp{i}, numNeigh(i)] = nstCollObj.getNeighbors(neuronNum(i));
               end
               if(min(numNeigh)==max(numNeigh)) %if all the same dimension
                   n=[];
                   for i=1:length(neuronNum)
                       n = [n;nTemp{i}];   %convert to a Matrix
                   end
               else %else return as a cell of neighbors
                   n=nTemp;
                   numNeighbors=max(numNeigh);
               end
           end
        end
       
        
        function setMinTime(nstCollObj,minTime)
            % setMinTime(nstCollObj,minTime)
            % calls setMinTime on all nspikeTrains in the collection and 
            % updates the minTime of the collection
            if(nargin<2)
                minTime=nstCollObj.minTime;
            end
            
            
            for i=1:nstCollObj.numSpikeTrains
                nstCollObj.nstrain{i}.setMinTime(minTime);
            end  
            nstCollObj.minTime = minTime;
%             nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end        
        function setMaxTime(nstCollObj,maxTime)
            % setMaxTime(nstCollObj,maxTime)
            % calls setMaxTime on all nspikeTrains in the collection and 
            % updates the maxTime of the collection
            if(nargin<2)
                maxTime=nstCollObj.maxTime;
            end
            
            for i =1:nstCollObj.numSpikeTrains
                nstCollObj.nstrain{i}.setMaxTime(maxTime);
            end
            nstCollObj.maxTime = maxTime;
%             nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end
        function setMask(nstCollObj,mask)
            if(length(mask)==nstCollObj.numSpikeTrains)
                if(max(mask)>1) %then these are indices
                    newMask = ones(1,length(mask));
                    nstCollObj.setNeuronMask(newMask);
                else %they are selectors of 1's and 0's
                    nstCollObj.setNeuronMask(mask);
                end
            else
                nstCollObj.setNeuronMaskFromInd(mask)
            end
        end
        function setNeuronMaskFromInd(nstCollObj,mask)
            % setNeuronMaskFromInd(nstCollObj,mask)
            % mask is a vector of indices into the nstColl.
            % all the nspikeTrains corresponding to the indices present in the mask 
            % will be visible in the collection
            newMask = zeros(1,nstCollObj.numSpikeTrains);
            newMask(mask) = 1;
            nstCollObj.setNeuronMask(newMask);
        end
        
        function setNeuronMask(nstCollObj,mask)
            % setNeuronMask(nstCollObj,mask)
            % mask needs to vector with binary entries
            % entries marked by 1s are visible, 0's are not visible
            if(length(mask)==nstCollObj.numSpikeTrains)
                nstCollObj.neuronMask = mask;
            end
        end
        function setNeighbors(nstCollObj,neighborArray)
            if(nargin<2)
                for i =1:nstCollObj.numSpikeTrains
                    neuronList = 1:nstCollObj.numSpikeTrains;
                    neighborArray(i,:)=find(i~=neuronList); %all but ith neuron
                end
            end
            if(~isempty(neighborArray))
                [numRows, ~]=size(neighborArray);
                if(numRows==nstCollObj.numSpikeTrains)
                    nstCollObj.neighbors = neighborArray;
                else
                    display('Neighbor Array is not of appropriate dimensions');
                end
            end
        end
        
        function ind=getIndFromMask(nstCollObj)
            % ind=getIndFromMask(nstCollObj)
            % ind is a row vector containing the indices of the currently
            % visible nspikeTrains in the collection
            
            ind=find(nstCollObj.neuronMask==1);
        end
        
        function ind=getIndFromMaskMinusOne(nstCollObj,neuron)
           ind = nstCollObj.getIndFromMask;
           ind = ind(ind~=neuron);
        end
        function answer = isNeuronMaskSet(nstCollObj)
            % answer = isNeuronMaskSet(nstCollObj)
            % answer =1 if any element of the neuronMask is set to zero.
            % answer =0 if all elements are set to 1.
            if(any(nstCollObj.neuronMask==0))
               answer =1;
            else 
               answer=0;
            end
        end        
         function answer = areNeighborsSet(nstCollObj)
            answer=~isempty(nstCollObj.neighbors);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function restoreToOriginal(nstCollObj,rMask)
            % restoreToOriginal(nstCollObj,rMask)
            % restores all elements of the nstColl to their original
            % states. If rMask is not specified then the neuronMask is
            % retained. If rMask =1, then the neuronMask is reset as well.
            if(nargin<2)
                rMask=0;
            end
            minTime=inf;
            maxTime=-inf;
            for i=1:nstCollObj.numSpikeTrains
                nstCollObj.getNST(i).restoreToOriginal;
                minTime=min(nstCollObj.getNST(i).minTime,minTime);
                maxTime=max(nstCollObj.getNST(i).maxTime,maxTime);
            end
            nstCollObj.setMinTime(minTime);
            nstCollObj.setMaxTime(maxTime);
            nstCollObj.sampleRate = nstCollObj.findMaxSampleRate;
            nstCollObj.enforceSampleRate;
            if(rMask==1)
                nstCollObj.resetMask;
            end 
%             nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end    
        function mSR= findMaxSampleRate(nstCollObj)
            mSR=-inf;
            for i=1:nstCollObj.numSpikeTrains
                mSR=max(nstCollObj.getNST(i).sampleRate,mSR);
            end
        end
        function resetMask(nstCollObj)
            %resetMask(nstCollObj)
            % resets the neuronMask to all 1's (All data visible)
            nstCollObj.neuronMask=ones(1,nstCollObj.numSpikeTrains);
%             nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end        
        function addToColl(nstCollObj,nst)
            % addToColl(nstCollObj,nst)
            % adds the nspikeTrain nst to the collection nstCollObj
            % nst: can be a single nspikeTrain or a cell array of 
            %      nspikeTrains
              if(isa(nst,'cell'))
                numElements=length(nst);
                for i=1:numElements
                    if(isa(nst{i},'nspikeTrain'))
                        nstCollObj.addSingleSpikeToColl(nst{i})
                    else
                        error('nstColl requires a cell array of nspikeTrain class elements');
                    end
                end
              elseif(isa(nst,'nspikeTrain'))
                  nstCollObj.addSingleSpikeToColl(nst);
              else
                  error('Can only add single spikes or cells with spikes');
              end
%               nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end        

        
        function uniqueNames = getUniqueNSTnames(nstCollObj,selectorArray)
            if(nargin<2)
                selectorArray = find(nstCollObj.neuronMask);
            end
            uniqueNames = unique(nstCollObj.neuronNames(selectorArray));
        end
        function names = getNSTnames(nstCollObj,selectorArray)
            if(nargin<2)
                selectorArray = find(nstCollObj.neuronMask);
            end
            names=nstCollObj.neuronNames(selectorArray);
            
            
        end
        
        
        function indices = getNSTIndicesFromName(nstCollObj,neuronName)
            if(nargin<2)
                neuronName = nstCollObj.getUniqueNSTnames;
            end
                
            if(isa(neuronName,'char'))
                names = nstCollObj.getNSTnames;
                indices = find(strcmp(names,neuronName));
            elseif(isa(neuronName,'cell'))
               for i=1:length(neuronName)
                  indices{i} = nstCollObj.getNSTIndicesFromName(neuronName{i});
               end
            end
            
        end
        
        function name = getNSTnameFromInd(nstCollObj,ind)
           if(ind>0 && nstCollObj.numSpikeTrains)
               name = nstCollObj.neuronNames(ind);
           else
               error('Index is out of bounds!');
           end
           
        end
        function nst = getNSTFromName(nstCollObj,neuronName)
            if(nargin<2)
                neuronName = nstCollObj.getUniqueNSTnames;
            end
                
            indices =nstCollObj.getNSTIndicesFromName(neuronName);
            nst = nstCollObj.getNST(indices); 
        end
        
        function nst = getNST(nstCollObj,index)
        % nst = getNST(nstCollObj,index)
        % if a nspikeTrain exists in the collection with the specified
        % index, it is returned in nst. 
        % If length(index)>1 then get a cell-array of nspikeTrains.
        % Care should be taken because it is
        % the actual nst and not a copy (recall handle behavior of
        % nspikeTrains)
            if(all(index>0) && all(index<=nstCollObj.numSpikeTrains))
                if(length(index)>1)
                    for i=1:length(index)
                        sampleRate = nstCollObj.sampleRate;
                        nst{i} = nstCollObj.nstrain{index(i)};%.nstCopy;
                        if(nst{i}.sampleRate~=sampleRate)
                            nst{i} = nst{i}.resample(sampleRate);
                        end
                    end
                else
                    sampleRate = nstCollObj.sampleRate;
                    nst = nstCollObj.nstrain{index};%.nstCopy;
                    if(nst.sampleRate~=sampleRate)
                        nst = nst.resample(sampleRate);
                    end
                end
            else
                error(['Index', num2str(index), ' out of bounds']);
            end
        end
        function resample(nstCollObj,sampleRate)
            %resample(nstCollObj,sampleRate)
            %resamples all the nspikeTrains in the collection at the
            %speciied sampleRate.
            %Updates teh sampleRate of the collection to the specified
            %value.
            
            if(nstCollObj.sampleRate~=sampleRate)
                for i=1:nstCollObj.numSpikeTrains
                    tempNST = nstCollObj.getNST(i);
                    tempNST.resample(sampleRate);
                    tempNST.setMinTime(nstCollObj.minTime);
                    tempNST.setMaxTime(nstCollObj.maxTime);

                end
                nstCollObj.sampleRate=sampleRate;
            end
%             nstCollObj.isSigRepBinary = nstCollObj.BinarySigRep;
        end
        function answer=isSigRepBinary(nstCollObj)
           answer = nstCollObj.BinarySigRep; 
        end
        function answer=BinarySigRep(nstCollObj)
            % answer=BinarySigRep(nstCollObj)
            % answer=1 if the all of the SignalObj representations of each
            % of the nspikeTrains in the collection have a binary
            % representation with the current collection parameters.
            % answer=0 if at least one nspikeTrain does not.
            tempAns = zeros(1,nstCollObj.numSpikeTrains);
            for i =1:nstCollObj.numSpikeTrains
                tempAns(i) = nstCollObj.getNST(i).isSigRepBinary;
            end
            
            answer=all(tempAns);

        end
        
        function ensembleCovariates = getEnsembleNeuronCovariates(nstCollObj,neuronNum,neighborIndex,windowTimes)
            % ensembleCovariates = getEnsembleNeuronCovariates(nstCollObj,neuronNum,neighborIndex,windowTimes)
            % returns a collection of covariates. Each covariate has number
            % of dimensions same as the number of history windows. There
            % will be the same number of covariates as neighbors.
     
            if(nargin<4 || isempty(windowTimes))
               windowTimes = [0 0.001];
            end
            
            if(nargin<3 || isempty(neighborIndex))
                allNeighbors=nstCollObj.getNeighbors(neuronNum);
            else
                allNeighbors=neighborIndex;
            end
            
            if(isa(windowTimes,'History'))
                histObj = windowTimes; %we were passed an actual history object;
            elseif(isa(windowTimes,'double'))
                histObj = History(windowTimes); % we got windowTimes;
            end
               
            ensembleCovariates = histObj.computeHistory(nstCollObj.getNST(1:nstCollObj.numSpikeTrains));
            ensembleCovariates.maskAwayAllExcept(allNeighbors);
            
            nstCollObj.addNeuronNamesToEnsCovColl(ensembleCovariates);
            
        end
        
        function addNeuronNamesToEnsCovColl(nstCollObj,ensembleCovariates)
            % addNeuronNamesToEnsCovColl(nstCollObj,ensembleCovariates)
            % Given a covariate collection of the ensemble effects, adds
            % the neuron number to each data dimension of each ensembleCovariate element. 
            
           
            for i=1:ensembleCovariates.numCov
                tempCov = ensembleCovariates.covArray{i};
                dataLabels = cell(1,tempCov.dimension);
                for j=1:tempCov.dimension;
                    name = nstCollObj.getNST(i).name;
                    if(isnumeric(name))
                        if(i>0 && i<10)
                            name = strcat(num2str(0),name);
                        end
                        dataLabels{j} = strcat('N',name, ':', tempCov.dataLabels{j});
                    else
                        dataLabels{j} = strcat(name, ':', tempCov.dataLabels{j});
                    end
                end
                tempCov.setDataLabels(dataLabels);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Change of Representation Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function dataMat = dataToMatrix(nstCollObj,selectorArray, binwidth,minTime,maxTime)
            % dataMat = dataToMatrix(nstCollObj,selectorArray, binwidth,minTime,maxTime)
            % dataMat is a matrix with number of rows corresponding to the
            % sampleRate of the nstColl and the minTime & maxTime values
            % that are currently set. The number of columns defaults to the
            % number of visible nspikeTrains as determined by the
            % neuronMask. This default behavior can be changed by specifyin
            % values for: selectorArray, binwidth, minTime, and maxTime
            % when calling the function
            
            if(nargin<5)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<4)
                minTime=nstCollObj.minTime;
            end
            if(nargin<3)
                binwidth=1/nstCollObj.sampleRate;
            end
            if(nargin<2)
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
             
            if(isa(selectorArray,'double'))
                if(any(selectorArray > nstCollObj.numSpikeTrains) || any(selectorArray < 1))
                    error('Neuron index is out of bounds');
                end
            elseif(isa(selectorArray,'cell')||isa(selectorArray,'char'));
               selectorArray = nstCollObj.getNSTIndicesFromName(selectorArray);
               if(isa(selectorArray,'cell'))
                   tempSelect = [];
                   for i=1:length(selectorArray)
                      tempSelect = [tempSelect selectorArray{i}];
                      
                   end
               end
            end
            
            %minTime and maxTime can be any values. We return a 
            %a matrix that starts at t=minTime and ends at t=maxTime;
%             maxTime
%             minTime
%             binwidth
            dataMat=zeros(floor(abs(maxTime-minTime)/binwidth)+1,length(selectorArray));
            testSig = nstCollObj.getNST(1).getSigRep(binwidth,minTime, maxTime);
            dataMat=zeros(length(testSig.dataToMatrix),length(selectorArray));
            for i=1:length(selectorArray) 
                nspikeSigRep=nstCollObj.getNST(selectorArray(i)).getSigRep(binwidth,minTime, maxTime);
%                   size(nspikeSigRep.dataToMatrix)
%                   size(dataMat)
                %A=nspikeSigRep.dataToMatrix;
                %dataMat(:,i)=A(1:length(dataMat),:);
%                 sum(nspikeSigRep.data)==size(nstCollObj.getNST(selectorArray(i)).spikeTimes,2)
                dataMat(:,i)=nspikeSigRep.dataToMatrix;
                d(i)=length(nstCollObj.getNST(selectorArray(i)).spikeTimes)-sum(dataMat(:,i));
            end
        end
        function spikeTrainObj = toSpikeTrain(nstCollObj,selectorArray,minTime,maxTime,windowTimes)
            if(nargin<5)
                windowTimes=[];
            end
            if(nargin<4 || isempty(maxTime))
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<3 || isempty(minTime))
                minTime=nstCollObj.minTime;
            end
            if(nargin<2 || isempty(selectorArray))
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            
            
             spikeTimes=[];
             offset=0;
             delta = 1/nstCollObj.sampleRate;
            nst=cell(1,length(selectorArray));
            
            if(isempty(windowTimes))
                for i=1:length(selectorArray)
                   if(i==1)
                       nst{i}=  nstCollObj.getNST(selectorArray(i));
                       spikeTimes = nst{i}.spikeTimes;
                       name = nst{i}.name;
                   else
                       offset = offset+nst{i-1}.maxTime + delta;
                       nst{i}=  nstCollObj.getNST(selectorArray(i));
                       if(~isempty(nst{i}.spikeTimes))
                        spikeTimes = [spikeTimes, nst{i}.spikeTimes+offset];
                       end

                   end
                end
            else
                if(length(selectorArray)==length(windowTimes)-1)
                    for i=1:length(selectorArray)
                       minTime=windowTimes(i);
                       deltaTW=windowTimes(i+1)-minTime;
                       if(i==1)
                           nst{i}=  nstCollObj.getNST(selectorArray(i));
                           spikeTimes = nst{i}.spikeTimes*deltaTW+minTime;
                           name = nst{i}.name;
                       else
    %                        offset = offset+nst{i-1}.maxTime + delta;
                           nst{i}=  nstCollObj.getNST(selectorArray(i));
                           spikeTimes = [spikeTimes, nst{i}.spikeTimes*deltaTW+minTime];

                       end
                    end
                else
                    error('Window Times must be 1 row longer than selectorArray');
                end
            end
            
            spikeTrainObj = nspikeTrain(spikeTimes);
            spikeTrainObj.setName(name); 
            maxTimeTot = maxTime*length(selectorArray);
            spikeTrainObj.setMinTime(minTime);
            spikeTrainObj.setMaxTime(maxTimeTot);
            spikeTrainObj.resample(1/delta);
            
        end
        function psthSignal = psth(nstCollObj,binwidth,selectorArray, minTime,maxTime)
            % psthSignal = psth(nstCollObj,selectorArray, binwidth,minTime,maxTime)
            % Given a collection a neural spike trains, sums the activity
            % in each time bin across all the neurons to generate a time
            % histogram of the neural firing. 
            % selectoryArray: can be used specify a subset of all the
            %                 neurons in the collection
            % binwidth: the size of the time bins used to sum the neural
            %           activity. Default is 1ms if not specified.
            % minTime:  time to start the psth. Default is the minTime of
            %           the nstCollObj.
            % maxTime:  time to end the psth. Degault to maxTime of
            %           nstCollObj.
            if(nargin<5)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<4)
                minTime=nstCollObj.minTime;
            end
            if(nargin<3 || isempty(selectorArray))
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            if(nargin<2)
                binwidth=.100; %100ms
            end
            windowTimes = minTime:binwidth:maxTime;
            if(~any(windowTimes==maxTime))
                windowTimes=[windowTimes, maxTime];
            end
                
            psthData=zeros(1,length(windowTimes));
                %dataMat = nstCollObj.dataToMatrix(nstCollObj,selectorArray,binwidth,minTime,maxTime);
            for i=1:length(selectorArray)    
                spikeTimes = nstCollObj.getNST(selectorArray(i)).getSpikeTimes;
                if(~isempty(spikeTimes))
                    
                    psthData = psthData+histc(spikeTimes,windowTimes);
                end
            end
                %due to how histc works add the last observation to the
                %previous one;
                tempPsth = psthData(1:end-1);
%                 tempPsth(end) = tempPsth(end)+psthData(end);
                psthData = tempPsth;
%                 unitPulseBasis=nstCollObj.generateUnitImpulseBasis(binwidth,minTime,maxTime);
%                 psthData = (unitPulseBasis.data*psthData')./binwidth;
                psthData = psthData./binwidth./length(selectorArray);
%                 unitPulseBasis=nstCollObj.generateUnitImpulseBasis(binwidth,minTime,maxTime,nstCollObj.sampleRate);
                
%                 psthData = unitPulseBasis.data*psthData';
%                 time = unitPulseBasis.time;
                time = (windowTimes(2:end)+windowTimes(1:end-1))/2;
                    
                
                psthSignal = SignalObj(time, psthData, 'PSTH','time','s','Hz');
%                 psthSignal.setMinTime(minTime);
%                 psthSignal.setMaxTime(maxTime);
               
%                 figure; bar(time,psth); ylabel('Spikes'); xlabel('time [sec]');
        end

        function psthSignal = psthBars(nstCollObj,binwidth,selectorArray, minTime,maxTime)
            % psthSignal = psthBars(nstCollObj,selectorArray, binwidth,minTime,maxTime)
            % Given a collection a neural spike trains, sums the activity
            % in each time bin across all the neurons to generate a time
            % histogram of the neural firing. Uses the implementation Baysian Adaptive Regression
            % Splines (developed by Wallstrom, Leibner and Kass) available
            % on Ryan C. Kelly's website:
            % <a href="matlab:web('http://www.cnbc.cmu.edu/~rkelly/code.html', '-helpbrowser')">Ryan C. Kelly's BARS for Matlab</a> 
            % Requires that: barsP.m, defaultParams.m, and nlsd_mex.c be in
            % the Matlab path and that nlsd_mex.c be compiled to the
            % corresponding mex file.
            % 
            % 
            % selectoryArray: can be used specify a subset of all the
            %                 neurons in the collection
            % binwidth: the size of the time bins used to sum the neural
            %           activity. Default is 1ms if not specified.
            % minTime:  time to start the psth. Default is the minTime of
            %           the nstCollObj.
            % maxTime:  time to end the psth. Degault to maxTime of
            %           nstCollObj.
            if(nargin<5)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<4)
                minTime=nstCollObj.minTime;
            end
            if(nargin<3 || isempty(selectorArray))
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            if(nargin<2)
                binwidth=.100; %100ms
            end
            time = minTime:binwidth:maxTime;
            psthData=zeros(1,length(time));
                %dataMat = nstCollObj.dataToMatrix(nstCollObj,selectorArray,binwidth,minTime,maxTime);
            for i=1:length(selectorArray)    
                spikeTimes = nstCollObj.getNST(selectorArray(i)).getSpikeTimes;
                if(~isempty(spikeTimes))
                    psthData = psthData+histc(spikeTimes,time);
                end
            end
            psthData = psthData./binwidth;
            psthData = psthData/length(selectorArray);
                %due to how histc works add the last observation to the
                %previous one;

            bp = defaultParams;
            bp.prior_id = 'POISSON';
            bp.dparams = 4;
  
  
            numTrials = length(selectorArray);
            fit = barsP(psthData,[minTime maxTime],numTrials);
            psthSignal = SignalObj(time, [fit.mode fit.mean fit.confBands], 'PSTH_{bars}','time','s','Hz',{'mode','mean','ciLower','ciUpper'});

        end

        
        
        function [xK,WK, Qhat,gammahat,logll,fitResults] = ssglm(nstCollObj,windowTimes,numBasis,numVarEstIter,fitType)
            if(nargin<5 || isempty(fitType))
               fitType='poisson';%'binomial'; 
            end
            if(nargin<4 || isempty(numVarEstIter))
                numVarEstIter=10;    
            end
            if(nargin<3 || isempty(numBasis))
                basisWidth=.02;
                numBasis = (nstCollObj.maxTime-nstCollObj.minTime)./basisWidth;
                
            end
            if(nargin<2)
               windowTimes = []; 
            end



            
            



            dN=nstCollObj.dataToMatrix';
            dN(dN>1)=1;
            basisWidth=(nstCollObj.maxTime-nstCollObj.minTime)/numBasis;
            [~, ~, psthResult] =nstCollObj.psthGLM(basisWidth,windowTimes,fitType);
            gamma0=psthResult.getHistCoeffs';%+.1*randn(size(histCoeffs));
            gamma0(isnan(gamma0))=-5;
            x0=psthResult.getCoeffs;

            
            Q0 = nstCollObj.estimateVarianceAcrossTrials(numBasis,windowTimes,numVarEstIter,fitType);

            if(any(diag(Q0)==0))
                Q0=Q0+0.001*diag(rand(numBasis,1));
            end
            A=eye(numBasis,numBasis);
%             fitType='poisson';
            delta = 1/nstCollObj.sampleRate;
            [xK,WK, Qhat,gammahat,logll]=DecodingAlgortihms.PPSS_EM(A,Q0,x0,dN,fitType,delta,gamma0,windowTimes, numBasis);

             minTime=nstCollObj.minTime; maxTime = nstCollObj.maxTime;
             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstCollObj.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end

            nCopy = nstCollObj.toSpikeTrain;
             

            histObj = History(windowTimes);
            cnt=1; [K,N]=size(dN);
            R=size(xK,1);
%             R=numBasis;
            clear beta otherLabels;
            lambdaData=[];
            nst=nstCollObj.nstrain;
            for k=1:K
                Hk{k}=histObj.computeHistory(nst{k}).dataToMatrix';
                
                stimK=basisMat*xK(:,k);

                histEffect=exp(gammahat(end,:)*Hk{k})';
                stimEffect=exp(stimK);
                lambdaDelta = histEffect.*stimEffect;
                
                lambdaData = [lambdaData;lambdaDelta/delta];
                
                for r=1:R
                    otherLabels{cnt} = ['b_{' num2str(r) ',' num2str(k) '}'];
                    beta(cnt) = xK(r,k);
                    cnt=cnt+1;
                    
                end
            end
            lambdaTime = minTime:delta:(length(lambdaData)-1)*delta;
            nCopy.setMaxTime(max(lambdaTime));
            nCopy.setMinTime(min(lambdaTime));
%             otherLabels  = tObj.getLabelsFromMask(neuronNumber);
            numLabels = length(otherLabels);
            histLabels  = histObj.computeHistory(nst{1}).getCovLabelsFromMask;
            otherLabels((numLabels+1):(numLabels+length(histLabels)))=histLabels;
            clear labels distrib stats b XvalData XvalTime;
            labels{1}  = otherLabels; % Labels change depending on presence/absense of History or ensCovHist
            numHist = length(histObj.windowTimes)-1;
%             histObj = tObj.history;
            ensHistObj = [];
%             [lambdaTemp, bTemp, devTemp, statsTemp,AICTemp,BICTemp,distribTemp] = Analysis.GLMFit(tObj,neuronNumber,i,Algorithm);
           lambdaIndexStr=1;
                lambda=Covariate(lambdaTime,lambdaData,...
                       '\Lambda(t)','time',...
                       's','Hz',strcat('\lambda_{',lambdaIndexStr,'}'));
                AIC = 2*length(otherLabels)-2*logll(end);
                BIC = -2*logll(end)+length(otherLabels)*log(length(lambdaData));
                statsTemp{1}=[];
                dev=-2*logll(end);
% lambda{i} = lambdaTemp; 
            b{1} = [beta';gammahat(end,:)']; 
            stats{1} = statsTemp;
%             dev(i) = devTemp;  
%             AIC(i)= AICTemp; 
%             BIC(i)= BICTemp;

            distrib{1} =fitType;
            currSpikes=nst;%nspikeColl.getNST(tObj.getNeuronIndFromName(neuronNames));
            for n=1:length(currSpikes)
                currSpikes{n} = currSpikes{n}.nstCopy;
                currSpikes{n}.setName(nCopy.name);
            end
            XvalData{1} = [];
            XvalTime{1} = [];
            spikeTraining = currSpikes;

                        
            fitResults=FitResult(nCopy,labels,numHist,histObj,ensHistObj,lambda,b, dev, stats,AIC,BIC,configColl,XvalData,XvalTime,distrib);
            DTCorrection=1;
            makePlot=0;
            Analysis.KSPlot(fitResults,DTCorrection,makePlot);
            Analysis.plotInvGausTrans(fitResults,makePlot);
            Analysis.plotFitResidual(fitResults,[],makePlot);
                %fitResults.computePlotParams;
            
        
            
        end
        
        
        function [psth, histSignal, psthResult] = psthGLM(nstCollObj, basisWidth,history,fitType,selectorArray,minTime,maxTime,sampleRate)
            if(nargin<8)
                sampleRate = 1/basisWidth;
            end
            if(nargin<7 || isempty(maxTime))
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<6 || isempty(minTime))
                minTime=nstCollObj.minTime;
            end
            if(nargin<5 || isempty(selectorArray))
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            if(nargin<4 ||isempty(fitType))
               fitType = 'poisson'; %'binomial';
            end
                
            if(nargin<3)
                history =[];
            end
            
            if(nargin<2)
                basisWidth = .100; %100ms
                
            end
            numBasis=ceil((maxTime-minTime)/basisWidth); 
            alphaVal = .05;
            unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);

            cc = CovColl({unitPulseBasis});
            trial = Trial(nstCollObj,cc);
            trial.setMinTime(minTime);
            trial.setMaxTime(maxTime);
%             for i=1:trial.nspikeColl.numSpikeTrains
%                 if(i==1)
%                     Y=num
%                 
   
            clear c;
            selfHist = history ; NeighborHist = []; sampleRate = nstCollObj.sampleRate; 
            LabelSelect = cell(1,unitPulseBasis.dimension+1);
            LabelSelect{1} = unitPulseBasis.name;
            LabelSelect(2:end) = unitPulseBasis.dataLabels(1:end);
            tc{1} = TrialConfig({LabelSelect},sampleRate,selfHist,NeighborHist); 
            if(~isempty(selfHist))
                tc{1}.setName('GLM-PSTH+Hist');
            else
                tc{1}.setName('GLM-PSTH');
            end
            cfgColl= ConfigColl(tc);
            warning off;
            
            if(strcmp(fitType,'poisson'))
                Algorithm='GLM';
            else
                Algorithm='BNLRCG';
            end
            
            batchMode =1;
            psthResult = Analysis.RunAnalysisForAllNeurons(trial,cfgColl,0,Algorithm,[],batchMode);
%             cfgColl.setConfig(trial,1);
%             for i=1:trial.nspikeColl.numSpikeTrains
%                 if(i==1)
%                     Y=trial.getSpikeVector(i);
%                     X=trial.getDesignMatrix(i);
%                 else
%                     Y=[Y;trial.getSpikeVector(i)];
%                     X=[X;trial.getDesignMatrix(i)];
%                 end
%             end
%             
% %                       y=tObj.getSpikeVector(neuronNumber);
% %             X=tObj.getDesignMatrix(neuronNumber);
% 
%             
%              distribution = 'poisson';
%              linkfunction = 'log';
%             
%              [bVals,dev,stats] = glmfit(X,Y,distribution, 'link', linkfunction,'constant','off');

            bVals = psthResult.b{1};
            stats = psthResult.stats{1};
            histVals = bVals(unitPulseBasis.dimension+1:end);
            statsHist.se = stats.se(unitPulseBasis.dimension+1:end);
            bVals=bVals(1:unitPulseBasis.dimension);
            statsBVals.se = stats.se(1:unitPulseBasis.dimension);
                
            % Could do multinomial fit here but not at this time...
%                 X=unitPulseBasis.data;
%                 Y=nstCollObj.dataToMatrix;
%                 [B,dev,stats]=mnrfit(X,Y);
            
            
%             results = Analysis.RunAnalysisForAllNeurons(trial,cfgColl,0);
%             warning on;
%             bVals = zeros(unitPulseBasis.dimension,length(results));
%             for i=1:length(results)
%                bVals(:,i) = results{i}.b{1}; 
%             end
            if(strcmp(fitType,'poisson'))
                expbVals = exp(bVals);
                expHistVals=exp(histVals);
            else
                expbVals = exp(bVals)./(1+exp(bVals));
                expHistVals = exp(histVals)./(1+exp(histVals));
            end
            
            unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
            psthData=(unitPulseBasis.data*expbVals)*(nstCollObj.sampleRate); % in Hz
%             psthData=(expbVals)*(nstCollObj.sampleRate); % in Hz
            %Remove Outliers%
            psthData(or(isnan(psthData),isinf(psthData)))=0;
%             modeVals = repmat(mode(psthData,2),[1 size(psthData,2)]);
            
%             psthData(psthData>nstCollObj.sampleRate*nstCollObj.numSpikeTrains)=0;
%           psthData=psthData*(basisWidth); %counts
            windowTimes = minTime:basisWidth:maxTime;
            
            if(~any(windowTimes==maxTime))
                windowTimes=[windowTimes, maxTime];
            end
            
            time = (windowTimes(2:end)+windowTimes(1:end-1))/2;
                       
            psth = Covariate(unitPulseBasis.time, psthData, 'PSTH_{glm}','time','s','Hz');
            
            
            
            mu=bVals;
            s=statsBVals.se;
            Mc=1000;
            for c=1:Mc
                z=normrnd(0,1,length(s),1);
                xKDraw(:,c)=mu+(s.*z);
            end
            if(strcmp(fitType,'poisson'))
%                 lower=logninv(alphaVal2/2,mu,s);
%                 lower=(unitPulseBasis.data*lower)*(nstCollObj.sampleRate); % in Hz
%                 upper=logninv(alphaVal2/2,mu,s);
%                 upper=(unitPulseBasis.data*upper)*(nstCollObj.sampleRate); % in Hz
            
                lambdaDraw=exp(xKDraw)*(nstCollObj.sampleRate);
                
                
            else %need to implement for binomial case
                lambdaDraw=exp(xKDraw)./(1+exp(xKDraw))*(nstCollObj.sampleRate);
                
            end
            lambdaDraw(isinf(lambdaDraw))=0;
            for k=1:length(s)
              [f,x] = ecdf(squeeze(lambdaDraw(k,:)));
              CIs(k,1) = x(find(f<alphaVal/2,1,'last'));
              CIs(k,2) = x(find(f>(1-alphaVal/2),1,'first'));
            end
           
%                 lower=logninv(0.025,mu,s);
            lower=(unitPulseBasis.data*CIs(:,1)); % in Hz
%                 upper=logninv(0.975,mu,s);
            upper=(unitPulseBasis.data*CIs(:,2)); % in Hz
  
            ciPSTHGLM = ConfidenceInterval(unitPulseBasis.time,[lower,upper],'CI_{psth_GLM}',psth.xlabelval,psth.xunits,psth.yunits);
            psth.setConfInterval(ciPSTHGLM);
            
            
            
            
          
            histTime=0:0.001:max(selfHist);
            if(~isempty(histTime))
                for i=1:length(selfHist)-1;
                    if(i==(length(selfHist)-1))
                        col = and(histTime>=selfHist(i),histTime<=selfHist(i+1))';
                    else
                        col = and(histTime>=selfHist(i),histTime<selfHist(i+1))';
                    end
                    basisMat(:,i) = col;
                end

    %             histVals = histVals(unitPulseBasis.dimension+1:end);

                histSignal = Covariate(histTime, basisMat*expHistVals, 'PSTH_{glm}','time','s','Hz');

                muH=histVals;
                sH=statsHist.se;
                clear xKDraw;
                Mc=1000;
                for c=1:Mc
                    z=normrnd(0,1,length(sH),1);
                    xKDraw(:,c)=sH.*z;
                end
            
                if(strcmp(fitType,'poisson'))
                   histDraw=exp(xKDraw)*(nstCollObj.sampleRate);
                
                
                else %need to implement for binomial case
                    histDraw=exp(xKDraw)./(1+exp(xKDraw))*(nstCollObj.sampleRate);
                
                end
                
                if(strcmp(fitType,'poisson'))
                    lowerH=logninv(0.025,muH,sH);
                    lowerH=(basisMat*lowerH); % in Hz
                    upperH=logninv(0.975,muH,sH);
                    upperH=(basisMat*upperH); % in Hz
                else %need to implement in binomial case
                    lowerH=logninv(0.025,muH,sH);
                    lowerH=(basisMat*lowerH); % in Hz
                    upperH=logninv(0.975,muH,sH);
                    upperH=(basisMat*upperH); % in Hz
                end
                
                for k=1:length(sH)
                    [f,x] = ecdf(squeeze(histDraw(k,:)));
                    CIsH(k,1) = x(find(f<alphaVal/2,1,'last'));
                    CIsH(k,2) = x(find(f>(1-alphaVal/2),1,'first'));
                end
                lowerH=(basisMat*CIsH(:,1)); % in Hz
                upperH=(basisMat*CIsH(:,2)); % in Hz
                
                ciPSTHGLMHist = ConfidenceInterval(histTime,[lowerH,upperH],'CI_{psth_GLMHIST}',psth.xlabelval,psth.xunits,psth.yunits);
                histSignal.setConfInterval(ciPSTHGLMHist);
            else
                histSignal = [];
            end
            
%             psth.setMinTime(minTime);
%             psth.setMaxTime(maxTime);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot(nstCollObj,selectorArray, minTime, maxTime,handle)
            if(nargin<5)
                handle = gca;
            end
            if(nargin<4)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<3)
                minTime=nstCollObj.minTime;
            end
            if(nargin<2)
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            
            %plotHandle = figure(handle);
            dHeight = 1;
            yOffset = 0:1:length(selectorArray)-1;
            yOffset = yOffset + dHeight/2;

            for i=1:length(selectorArray)
                currentObj = nstCollObj.getNST(selectorArray(i));
                currentObj.plot(dHeight, yOffset(i), handle); hold all;
                yticklabels{i} = currentObj.name;
                if(strcmp(yticklabels{i},''));
                    yticklabels{i} = num2str(selectorArray(i));
                end
            end
            xlabel('$$ time [s] $$','Interpreter','latex');            
            set(gca,'YTick',yOffset);
            set(gca,'YTickLabel',yticklabels);
            axis tight;
            v=axis;
%             minTime
%             maxTime
%             
            axis([minTime,maxTime,v(3),v(4)]);

%             plotHandle = figure(handle);
%             binwidth = max(nstCollObj.getMaxBinSizeBinary,.001);
%             nMat=nstCollObj.dataToMatrix(selectorArray,binwidth,minTime,maxTime);
% %             time = nstCollObj.minTime:(1/nstCollObj.sampleRate):nstCollObj.maxTime;
%             nMat(nMat>1)=1;
%             nMat = ~nMat;
%             imagesc(flipud(nMat'));
%             colormap(gray);
%             xt=get(gca,'xtick')*binwidth;
%             xtStr= num2str(xt');
%             set(gca,'xtickLabel',xtStr);
%             
%             for i=1:length(selectorArray)
%                currentObj = nstCollObj.getNST(selectorArray(i));
%                yticklabels{i} = currentObj.name;
%                if(strcmp(yticklabels{i},''));
%                     yticklabels{i} = num2str(selectorArray(i));
%                end
%             end
%             set(gca,'YTick',selectorArray,'YTickLabel',fliplr(yticklabels));
%             xlabel('$$ time [s] $$','Interpreter','latex');       


        end
   
        function plotISIHistogram(nstCollObj,selectorArray, minTime, maxTime,handle)
            if(nargin<5)
                handle = gca;
            end
            if(nargin<4)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<3)
                minTime=nstCollObj.minTime;
            end
            if(nargin<2)
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            
                        
            for i=1:length(selectorArray)
                currentObj = nstCollObj.getNST(selectorArray(i));
                figure;
                currentObj.plotISIHistogram(minTime,maxTime); %hold all;
                yticklabels{i} = currentObj.name;
                if(strcmp(yticklabels{i},''));
                    yticklabels{i} = num2str(selectorArray(i));
                end
            end
            %xlabel('$$ time [s] $$','Interpreter','latex');            
            %set(gca,'YTick',yOffset);
            %set(gca,'YTickLabel',yticklabels);
            %axis tight;
            %v=axis;
%             minTime
%             maxTime
%             
            %axis([minTime,maxTime,v(3),v(4)]);
        end
        
        function plotExponentialFit(nstCollObj,selectorArray, minTime, maxTime,numBins,handle)
            if(nargin<6)
                handle = gca;
            end
            if(nargin<5 || isempty(numBins))
                numBins =[];
            end
            if(nargin<4)
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<3)
                minTime=nstCollObj.minTime;
            end
            if(nargin<2)
                if(nstCollObj.isNeuronMaskSet)
                    selectorArray = nstCollObj.getIndFromMask;
                else
                    selectorArray = 1:nstCollObj.numSpikeTrains;
                end
            end
            
                        
            for i=1:length(selectorArray)
                currentObj = nstCollObj.getNST(selectorArray(i));
                
                currentObj.plotExponentialFit(minTime,maxTime,numBins); %hold all;
                yticklabels{i} = currentObj.name;
                if(strcmp(yticklabels{i},''));
                    yticklabels{i} = num2str(selectorArray(i));
                end
            end
            %xlabel('$$ time [s] $$','Interpreter','latex');            
            %set(gca,'YTick',yOffset);
            %set(gca,'YTickLabel',yticklabels);
            %axis tight;
            %v=axis;
%             minTime
%             maxTime
%             
            %axis([minTime,maxTime,v(3),v(4)]);
        end
        function varEst = estimateVarianceAcrossTrials(nstCollObj,numBasis,windowTimes,numIter,fitType)
            %returns a estimate of the variance within each basis
            if(nargin<5 || isempty(fitType))
               fitType = 'poisson'; 
            end
            if(nargin<4 || isempty(numIter))
                numIter=20;
            end
            if(nargin<3 || isempty(windowTimes))
                windowTimes = [];
            end
            if(nargin<2 ||isempty(numBasis))
                numBasis = 20;
            end

            echo off;
            coeffs = zeros(numBasis,numIter);
            numRealizations = nstCollObj.numSpikeTrains;
            basisWidth = (nstCollObj.maxTime - nstCollObj.minTime)./numBasis;
            sumNumber=numRealizations/2-1;
            delta = 1/nstCollObj.sampleRate;
            minTime = nstCollObj.minTime;
            maxTime = nstCollObj.maxTime;
            

            for i=1:min([numIter/2 sumNumber])

                spikeCollTemp=nstColl(nstCollObj.getNST(i:i+sumNumber));


                spikeCollTemp.resample(1/delta);
                spikeCollTemp.setMaxTime(maxTime);
                spikeCollTemp.setMinTime(minTime);
                [~, ~, psthResultT] =spikeCollTemp.psthGLM(basisWidth,windowTimes,fitType);
                coeffs(:,i)=psthResultT.getCoeffs;
            end
            for i=numRealizations:-1:(numRealizations-min([numIter/2 sumNumber])+1)
                spikeCollTemp=nstColl(nstCollObj.getNST(i:-1:i-sumNumber));


                spikeCollTemp.resample(1/delta);
                spikeCollTemp.setMaxTime(maxTime);
                spikeCollTemp.setMinTime(minTime);
                [~, ~, psthResultT] =spikeCollTemp.psthGLM(basisWidth,windowTimes,fitType);
                coeffs(:,i)=psthResultT.getCoeffs;
            end 
            
            %Remove zero columns
            for i=1:size(coeffs,1)
                CoeffsTemp(i,:) = coeffs(i,coeffs(i,:)~=0);
            end
            
            coeffs=CoeffsTemp;
            NTerms=4;A=1; B=ones(1,NTerms)./NTerms;
            coeffs(isnan(coeffs))=0;
%             coeffs(exp(coeffs)/delta>10)=log(10*delta); %dont allow for really large coeffs
            if(size(coeffs',1)>3*NTerms)
                fcoeffs = filtfilt(B,A,coeffs')';
            else
                fcoeffs = coeffs;
            end
            
            varEst=nanvar(diff(fcoeffs,[],2),[],2);

%             varEst(varEst>.001)=0.0001; %avoid large estimates of the sample variance

            varEst=diag(varEst);
            echo on;
        end
        function windowedSpikeTimes=getSpikeTimes(nstCollObj, minTime, maxTime)
            if(nargin<3 || isempty(maxTime))
                maxTime=nstCollObj.maxTime;
            end
            if(nargin<2 || isempty(minTime))
                minTime=nstCollObj.minTime;
            end
            
            ind = nstCollObj.getIndFromMask;
            windowedSpikeTimes=cell(length(ind),1);
            for i=ind
                if(i==1)
                    count=1;
                end
                currSpike = nstCollObj.getNST(i);
                windowedSpikeTimes{count} = currSpike.getSpikeTimes;
                count=count+1;
            end
        end
        
        function structure = toStructure(nstCollObj)
            fnames = fieldnames(nstCollObj);
            nstCollObj.resetMask; %otherwise masked data will not get saved!!
            for i=1:length(fnames)
                currObj = nstCollObj.(fnames{i});
                if(isa(currObj,'double')||isa(currObj,'cell')||isa(currObj,'logical'))
                    if(strcmp(fnames{i},'nstrain'))
                        for j=1:nstCollObj.numSpikeTrains
                            structure.(fnames{i}){j} = nstCollObj.(fnames{i}){j}.toStructure;
                        end
                    else
                        structure.(fnames{i}) =  currObj;
                    end
                end
            end
        end
    end
    methods (Static)
        function nstCollObj = fromStructure(structure)
            nst = cell(1,structure.numSpikeTrains);
            for i=1:structure.numSpikeTrains;
                nst{i} = nspikeTrain.fromStructure(structure.nstrain{i});
            end
            nstCollObj = nstColl(nst);
            nstCollObj.setMinTime(structure.minTime);
            nstCollObj.setMaxTime(structure.maxTime);
            nstCollObj.setNeighbors(structure.neighbors);
        end
        function unitPulseBasis=generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate)
            % Samplerate determines number of samples per second for the
            % unit pulse functions. If basisWidth is larger than the sample
            % rate of nstCollObj, when a trial is formed to obtain the GLM
            % psth, the up-sampled pulses will no longer be square. Better
            % to upsample the spikeTrain instead.
            if(nargin<5)
                sampleRate=1000;
            end
             
            windowTimes = minTime:basisWidth:(maxTime);
            
            
            if(~any(windowTimes==maxTime))
                windowTimes=[windowTimes, maxTime];
            end
            numBasis=length(windowTimes)-1;   
%             timeVec = (minTime:(1/nstCollObj.sampleRate):maxTime)';
            timeVec = (minTime:(1/sampleRate):maxTime)';
            dataMat = zeros(length(timeVec),length(windowTimes)-1);
            
            for i=1:length(windowTimes)-1
               dataMat(:,i) = and(timeVec>=(windowTimes(i)),timeVec<windowTimes(i+1));
               if(i==length(windowTimes)-1)
                  dataMat(:,i) = and(timeVec>=(windowTimes(i)),timeVec<=windowTimes(i+1)); 
               else
                  dataMat(:,i) = and(timeVec>=(windowTimes(i)),timeVec<windowTimes(i+1));
               end
            end
            
            for i=1:numBasis
                if(i<10)
                    dataLabels{i} = strcat('b0',num2str(i));
                else
                    dataLabels{i} = strcat('b',num2str(i));
                end
            end
                
            unitPulseBasis = Covariate(timeVec,dataMat,'UnitPulseBasis','time','s','',dataLabels);
            
            
        end
        
        
    end
    methods (Access = private)
        function addSingleSpikeToColl(nstCollObj,nst)
          if(isa(nst,'nspikeTrain'))
                if(isempty(nst.name))
                    nst.setName(num2str(nstCollObj.numSpikeTrains+1));
                end
                nstCollObj.nstrain{nstCollObj.numSpikeTrains+1}= nst;
                nstCollObj.neuronNames{nstCollObj.numSpikeTrains+1} = nst.name;
                nstCollObj.updateTimes(nst);
                nstCollObj.numSpikeTrains = nstCollObj.numSpikeTrains + 1;
                nstCollObj.neuronMask = [nstCollObj.neuronMask,1];
                nstCollObj.sampleRate = max(nstCollObj.sampleRate,nst.sampleRate);
                nstCollObj.enforceSampleRate;

          else 
                error('Can only add neural spike trains to collection');
          end
        end  
        function ensureConsistancy(nstCollObj)
                nstCollObj.enforceSampleRate;
                nstCollObj.setMinTime;
                nstCollObj.setMaxTime;
        end        
        function enforceSampleRate(nstCollObj) 
            for i=1:nstCollObj.numSpikeTrains;
                currSpike = nstCollObj.getNST(i);
                if(currSpike.sampleRate~=nstCollObj.sampleRate)
                   currSpike.resample(nstCollObj.sampleRate);
                end
            end
        end       
        function updateTimes(nstCollObj,nst)
            %updateTimes(nstCollObj,nst)
            % Given the nspikeTrain nst, the times in the nstColl are
            % updated if nst.minTime is smaller than the collection minTime
            % or if nst.maxTime is larger than the collection maxTime
            if(nst.minTime<=nstCollObj.minTime)
                nstCollObj.setMinTime(nst.minTime);
            else
                nst.setMinTime(nstCollObj.minTime);
            end
            if(nst.maxTime>=nstCollObj.maxTime)
                nstCollObj.setMaxTime(nst.maxTime);
            else
                nst.setMaxTime(nstCollObj.maxTime);
            end
        end
        
      
        
    end
    
end

