classdef CovColl <handle 
% COVCOLL a collection of covariates. Allows multiple covariates that were
% recorded simultaneously to be treated as single unit. Operations such
% as resampling, setting time windows, etc can then be performed on the
% collection as a whole. 
%
% When covariates are accessed through the getCov function, copies of the covariates are return
% the original covariates remain intact. The covariate collection
% remembers the masked states, shifts, etc. so that these are applied
% to the signal right before it is returned.
% 
% <a href="matlab: methods('CovColl')">methods</a>
% <a href="matlab:web('CovCollExamples.html', '-helpbrowser')">CovColl Examples</a> 
%
% see also <a href="matlab:help('SignalObj')">SignalObj</a>, <a href="matlab:help('Covariate')">Covariate</a>
%
% Reference page in Help browser
% <a href="matlab: doc('CovColl')">doc CovColl</a>
% 

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
        covArray; %An array of covariate objects
        covDimensions; % at each position has the number of dimensions of each covariate
        numCov; % a running count of how many covariates are in object
        minTime; % Time data occurs
        maxTime; % Time last data point occurs in object
        covMask; % covariates that are currently selected
        covShift;% time lag for covariates
        sampleRate % sampleRate for all covariates
   end

   properties (Hidden)
        originalSampleRate;
        originalMinTime;
        originalMaxTime;
       
   end
       
    
    methods
        function ccObj=CovColl(cov,varargin)
            % ccObj=CovColl(cov,varargin)
            % Creates a collection of covariates from a cell array of
            % objects of the class Covariate <a href="matlab:help('Covariate')">Covariate</a>
             if(nargin<1)
                 cov=[];
             end
             ccObj.numCov = 0;
             ccObj.minTime=inf;
             ccObj.maxTime=-inf;
             ccObj.originalSampleRate=[];
             ccObj.originalMinTime = [];
             ccObj.originalMaxTime = [];
             ccObj.covArray=[];
             ccObj.covDimensions=[];
             ccObj.covMask = [];
             ccObj.covShift = 0;
             ccObj.addToColl(cov);
             if(nargin>1)
                 for i=1:length(varargin)
                    ccObj.addToColl(varargin{i});
                 end
             end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setMinTime(ccObj,minTime)
            % setMinTime(ccObj,minTime)
            % sets the minimum time for all the covariates in the
            % collection to minTime
            if(nargin<2 || isempty(minTime))
                minTime=ccObj.findMinTime;
            end
            if(isempty(ccObj.originalMinTime))
                 ccObj.originalMinTime=ccObj.minTime;
            end
%             for i=1:ccObj.numCov
%                 tempC = ccObj.covArray{i};
%                 tempC.setMinTime(minTime);
%             end
            ccObj.minTime=minTime;
           end
        function setMaxTime(ccObj,maxTime)
            % setMaxTime(ccObj,maxTime)
            % sets the maximum time for all the covariates in the
            % collection to maxTime
            
            if(nargin<2 || isempty(maxTime))
                maxTime=ccObj.findMaxTime;
            end
            if(isempty(ccObj.originalMaxTime))
                 ccObj.originalMaxTime=ccObj.maxTime;
            end
%             for i=1:ccObj.numCov
%                 tempC =ccObj.covArray{i};
%                 tempC.setMaxTime(maxTime);
%             end
            ccObj.maxTime = maxTime;
        end            
        function setSampleRate(ccObj, sampleRate)
            % setSampleRate(ccObj, sampleRate)
            % resample all of the covariates to the specified sampleRate
             if(isempty(ccObj.originalSampleRate))
                 ccObj.originalSampleRate=ccObj.sampleRate;
             end
%              minTime = ccObj.minTime;
%              maxTime = ccObj.maxTime;
             
             ccObj.sampleRate = sampleRate;
             ccObj.enforceSampleRate;
%              ccObj.restrictToTimeWindow(minTime,maxTime);
             
        end           
        function setMask(ccObj,cellInput)
            % setMask(ccObj,cellInput)
            % specify which covariates are to be used
            
             selectorCell = ccObj.generateSelectorCell(cellInput);
             ccObj.setMasksFromSelector(selectorCell);
             for i=1:ccObj.numCov
                cov=ccObj.getCov(i);
                cov.setMask(ccObj.covMask{i})
             end   
        end        
        function dataMask=getCovDataMask(ccObj,identifier)
            % dataMask=getCovDataMask(ccObj,identifier)
            % returns the dataMask for the covariate specified by
            % indentifier
            cov=ccObj.covArray{identifier};
            dataMask = cov.dataMask;
        end     
        function answer=isCovMaskSet(ccObj)
            % answer=isCovMaskSet(ccObj)
            % returns 1 if any Covariate has any component that is masked
            % away, otherwise returns 0.
            answer =0;
            for i=1:ccObj.numCov
                if(any(ccObj.covMask{i}==0))
                    answer =1;
                    break;
                end
            end
        end                
        function n=nActCovar(ccObj)
                % n=nActCovar(ccObj)
                % Returns the effective number of a covariates. Any
                % covariate with at least one unmasked component
                % contributes to n. Any covariate with all components
                % masked away is not counted.
                selectorArray = ccObj.getSelectorFromMasks;
                n=numActCov(selectorArray);
        end
            
        function maskAwayCov(ccObj,identifier)
            % maskAwayCov(ccObj,identifier) 
            % masks away all the components of the covariates specified by
            % indentifier
            cov=ccObj.getCov(identifier);
            if(isa(cov,'Covariate'))
                cov = {cov}; % make it a cell even if just one
            end
            for j=1:length(cov)
                covIndex = ccObj.getCovIndicesFromNames(cov{j}.name);            
                newMask = cell(1,ccObj.numCov);
                for i=1:ccObj.numCov
                    if(i==covIndex)
                        newMask{i} = zeros(1,length(ccObj.covMask{i}));
                    else
                        newMask{i} = ccObj.covMask{i};
                    end
                end
                ccObj.setMask(ccObj.getSelectorFromMasks(newMask));
            end
        end
        
        function ccObj2 = copy(ccObj)
            cov = cell(length(ccObj.numCov),1);
            for i=1:ccObj.numCov
               cov{i} = ccObj.getCov(i).copySignal;
            end
            ccObj2 = CovColl(cov);
            
            
        end
        function maskAwayOnlyCov(ccObj,identifier)
            % maskAwayOnlyCov(ccObj,identifier)
            % makes all components of all covariates visible and then masks
            % away the covariates specified by indentifier.
            ccObj.resetMask;
            ccObj.maskAwayCov(identifier);
        end
        
        function maskAwayAllExcept(ccObj, identifier)
            % maskAwayAllExcept(ccObj, identifier)
            % masks away all covariates except that specified by
            % identifier
            offset=0;
            maskList = zeros(1,ccObj.numCov - length(identifier));
            for i =1:ccObj.numCov
                if(~any(i==identifier))% i is not in any element of identifier
                    offset=offset+1;
                    maskList(offset) = i;
                end
            end
            ccObj.maskAwayOnlyCov(maskList);
        end
        
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function cov = getCov(ccObj, identifier)
            % cov = getCov(ccObj, identifier)
            % returns a single covariate if only one is requested.
            % Otherwise returns a cell array of covariates, one for each
            % identifier.
            % The identifier can be:
            %   doubles: specifying the number of the covariate in the
            %   collection
            %   strings: specifying the name of the covariate
            %   cell array of strings: specifying multiple covariates by
            %   their name.
            holdVals=1;
            if(isa(identifier,'double'))
                if(length(identifier)==1)
                    cov=ccObj.covArray{identifier}.copySignal;
                    cov.setMask(ccObj.covMask{identifier});
                    if(ccObj.covShift~=0)
                        cov=cov.shift(ccObj.covShift);
                    end
                    if(cov.minTime~=ccObj.minTime || cov.maxTime~=ccObj.maxTime)
                        cov=cov.getSigInTimeWindow(ccObj.minTime,ccObj.maxTime,holdVals);
                    end
%                    cov=cov.resample(ccObj.sampleRate);
                    
                else
                    cov=cell(1,length(identifier));
                    for i=1:length(identifier)
                        cov{i}=ccObj.getCov(identifier(i));
                        %cov{i}=ccObj.covArray{identifier(i)};
                    end
                end
            elseif(isa(identifier,'char'))
                %cov=ccObj.covArray{ccObj.getCovIndFromName(identifier)};
                 cov=ccObj.getCov(ccObj.getCovIndFromName(identifier));
            elseif(isa(identifier,'cell'))
                cov=cell(1,length(identifier));
                if(isa(identifier{1},'char'))
                    for i=1:length(identifier)
                        %cov{i}=ccObj.covArray{ccObj.getCovIndFromName(identifier{i})};
                         cov{i}=ccObj.getCov(identifier{i});
                    end
                else
                    error('Identifier cells must contain strings!');
                end
            end    
        end         
        function ind = getCovIndicesFromNames(ccObj,name)
            % ind = getCovIndicesFromNames(ccObj,name)
            % returns a vector of indices for each covariate name
            % specified.
            if(isa(name,'cell'))
                if(isa(name{1},'char'))
                    ind=zeros(1,length(name));
                    for i=1:length(name)
                        ind(i)=ccObj.getCovIndFromName(name{i});
                    end
                else
                    error('Cell must contain strings!');
                end
            elseif(isa(name,'char'))
                ind=ccObj.getCovIndFromName(name);
            else
                error('Need either cells with strings or a single string!');
            end
                       
        end 
        function dim = getCovDimension(ccObj,identifier)
            % dim = getCovDimension(ccObj,identifier)
            % returns a vector with the dimension of covariate i at
            % position i.
           covs = ccObj.getCov(identifier);
           dim = zeros(1,length(covs));
           for i=1:length(covs)
               dim(i)=covs{i}.dimension;
           end
        end      
        function l   = getAllCovLabels(ccObj)
            % l   = getAllCovLabels(ccObj)
            % returns a cell array of strings with the covariate names
            offset=0;
            l=cell(1,length(ccObj.flattenCovMask));
           for i=1:ccObj.numCov
               tempCov = ccObj.getCov(i);
               for j=1:tempCov.dimension
                   l{j+offset} = tempCov.dataLabels{j};
               end
               offset=offset+tempCov.dimension;
           end
        end
        
        function l = getCovLabelsFromMask(ccObj)
            % l = getCovLabelsFromMask(ccObj)
            % returns a list of all the the dataLabels that are currently
            % visible (i.e. unmasked).
            offset=0;
            l={};
            for i=1:ccObj.numCov
               tempCov = ccObj.getCov(i);
               for j=1:tempCov.dimension
                   if(ccObj.covMask{i}(j)==1)
                       offset=offset+1;
                        l{offset} = tempCov.dataLabels{j};
                   end
               end
            end    
            
        end
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function structure = toStructure(ccObj)
            fnames = fieldnames(ccObj);
            ccObj.resetMask; %otherwise masked data will not get saved!!
            for i=1:length(fnames)
                currObj = ccObj.(fnames{i});
                if(isa(currObj,'double')||isa(currObj,'cell'))
                    if(strcmp(fnames{i},'covArray'))
                        for j=1:ccObj.numCov
                            structure.(fnames{i}){j} = ccObj.(fnames{i}){j}.toStructure;
                        end
                    else
                        structure.(fnames{i}) =  currObj;
                    end
                end
            end
            
        end
        
        function minTime = findMinTime(ccObj)
            % minTime = findMinTime(ccObj)
            % finds the minimum minTime from all covariates
            minTime=inf;
            for i=1:ccObj.numCov
                minTime = min(ccObj.covArray{i}.minTime,minTime);
            end       
            minTime = minTime+ccObj.covShift;
        end
        function maxTime = findMaxTime(ccObj)
            % maxTime = findMaxTime(ccObj)
            % finds that maximum maxTime from all covariates.
            
            maxTime=-inf;
            for i=1:ccObj.numCov
                maxTime = max(ccObj.covArray{i}.maxTime+ccObj.covShift,maxTime);
            end   
            maxTime = maxTime+ccObj.covShift;
        end    

        function addToColl(ccObj,cov)
            % addToColl(ccObj,cov)
            % add one or several covariates to the current collection.
            % can specify cell of covariates, a single covariate, or a
            % covariate CovColl.
            if(~isempty(cov))
              if(isa(cov,'cell'))
                ccObj.addCovCellToColl(cov);
              elseif(isa(cov,'Covariate'))
                ccObj.addSingleCovToColl(cov);    
              elseif(isa(cov,'CovColl'));
                ccObj.addCovCollection(cov);
              else 
                error('Can only add covariates to CovColl');
              end
            end
            ccObj.enforceSampleRate;
        end          
        function addCovCollection(ccObj,cov)
            % addCovCollection(ccObj,cov)
            % adds a CovColl to the current collection
            covCell=cov.covArray;
            ccObj.addCovCellToColl(covCell);      
        end        

        
        function answer = isCovPresent(ccObj,cov)
            % answer = isCovPresent(ccObj,cov)
            % returns 1 if covariate is present in the CovColl.
            % inputs can be a covariate, a string corresponding to the name
            % of the covariate, or the number of the covariate in the
            % collection.
            if(isa(cov,'Covariate'))
                if(strcmp(cov.name,''))
                    display('Covariate does not have name');
                    answer=0;
                else
                    index=ccObj.getCovIndFromName(cov.name);
                    if(isempty(index))
                        answer = 0;
                    else
                        answer = 1;
                    end
                end
            elseif(isa(cov,'char'))
                covar=ccObj.getCov(cov);
                answer=ccObj.isCovPresent(covar);
            elseif(isa(cov,'double'))
                if((cov>0)&&(cov<ccObj.numCov))
                    answer=1;
                else
                    answer=0;
                end
            else
                error('Need either covariate class or name of covariate or index of covariate');
            end
        end
        
        function resample(ccObj,sampleRate)
            % resample(ccObj,sampleRate) 
            % resamples all the covariates in the collection to the new
            % sampleRate.
            ccObj.setSampleRate(sampleRate);
            ccObj.enforceSampleRate;
        end
        function restoreToOriginal(ccObj)
        % restoreToOriginal(ccObj)
        % returns the CovColl to the original minTime, maxTime, and
        % sampleRate. covShift is returned to zero.
%             minTime=inf;
%             maxTime=-inf;
              %minTime = ccObj.findMinTime;
              %maxTime = ccObj.findMaxTime;
%             for i=1:ccObj.numCov
%                 tempCov = ccObj.getCov(i);
%                 tempCov.restoreToOriginal;
%                 minTime=min(tempCov.minTime,minTime);
%                 maxTime=max(tempCov.maxTime,maxTime);
%             end
            ccObj.covShift = 0;
            ccObj.setSampleRate(ccObj.originalSampleRate);
            ccObj.setMinTime(ccObj.findMinTime);
            ccObj.setMaxTime(ccObj.findMaxTime);
            %ccObj.setMinTime(minTime);
            %ccObj.setMaxTime(maxTime);            
        end     
        function restrictToTimeWindow(ccObj,wMin,wMax)
            % restrictToTimeWindow(ccObj,wMin,wMax)
            % sets minTime to wMin, and maxTime to wMax
            ccObj.setMinTime(wMin);
            ccObj.setMaxTime(wMax);
%             for i=1:ccObj.numCov
%                 ccObj.getCov(i).setMinTime(wMin);
%                 ccObj.getCov(i).setMaxTime(wMax);
%             end
        end
        function removeCovariate(ccObj,identifier)
            % removeCovariate(ccObj,identifier)
            % removes the specified covariate from the collection
            ccObj.removeFromColl(identifier);
        end
        function resetMask(ccObj)
            % resetMask(ccObj)
            % makes all covariates visible
            for i=1:ccObj.numCov
                ccObj.covArray{i}.resetMask;
                ccObj.covMask{i}=ccObj.getCovDataMask(i);
            end           
        end        
        function enforceSampleRate(ccObj) 
            % enforceSampleRate(ccObj) 
            % makes sure that all covariates have the same sampleRate as
            % that in ccObj.sampleRate;
            for i=1:ccObj.numCov;
                currCov = ccObj.covArray{i}; %change the actual sample rate of the objects
                if(and(and(round(currCov.sampleRate*1000)/1000~=round(ccObj.sampleRate*1000)/1000,~isnan(currCov.sampleRate)),~isnan(ccObj.sampleRate)))
                   currCov.resampleMe(ccObj.sampleRate);
                end
            end
        end  
        
        function ccObj = setCovShift(ccObj, deltaT, identifier)
            % setCovShift(ccObj, deltaT, identifier)
            % Note: identifier currently not used
            % shifts ALL covariates by deltaT
            if(nargin<3)
                identifier=ccObj.getSelectorFromMasks;
            end
%             covars=ccObj.getCov(identifier);
%             for i=1:length(covars)
%                 covars{i}.shift(deltaT);
%             end
            ccObj.resetCovShift;
            ccObj.covShift=deltaT;
            ccObj.setMinTime(ccObj.minTime+deltaT); %make sure minTime is consistent
            ccObj.setMaxTime(ccObj.maxTime+deltaT); %make sure maxTime is consistent
           
        end
        
        function resetCovShift(ccObj)
            ccObj.covShift=0;
            ccObj.setMinTime; %make sure minTime is consistent
            ccObj.setMaxTime; %make sure maxTime is consistent
        end
        
        function flatMask = flattenCovMask(ccObj)
            covMask=ccObj.covMask;
            if(isa(covMask,'double'))
                flatMask=covMask;
            elseif(isa(covMask,'cell'))
                flatMask=[];
                for i=1:length(covMask)
                    flatMask = [flatMask covMask{i}];
                end
            else
                error('covMask must be either a cell or a double');
            end
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Change of Representation Functions 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dataMat = dataToMatrix(ccObj,repType,dataSelector,varargin)
            % dataMat = dataToMatrix(ccObj,repType,dataSelector,varargin)
            % returns the matrix representation of the CovColl. 
            % repType: 'standard' or 'zero-mean'
            % dataSelector: same as input to getCov
            if(nargin<3)
                dataSelector=ccObj.getSelectorFromMasks;
            end
            if(nargin<2)
                repType='standard';
            end
            
            if(ccObj.isaSelectorCell(dataSelector))
                dataMat=ccObj.dataToMatrixFromSel(repType,dataSelector,varargin{:});
            else %we assume these are names
                dataMat=ccObj.dataToMatrixFromNames(repType,dataSelector,varargin{:});
            end
        end            
        function dataMat = dataToMatrixFromNames(ccObj,repType,dataSelector,varargin)
            selectorCell=ccObj.generateSelectorCell(dataSelector);
            dataMat=ccObj.dataToMatrixFromSel(repType, selectorCell, varargin{:});
        end
        function dataMat = dataToMatrixFromSel(ccObj,repType, selectorCell,varargin)%, binwidth,minTime,maxTime)
%             if(nargin<6)
%                 maxTime=ccObj.maxTime;
%             end
%             if(nargin<5)
%                 minTime=ccObj.minTime;
%             end
%             if(nargin<4)
%                 binwidth=1/ccObj.sampleRate;
%             end
            if(nargin<3)
                if(ccObj.isCovMaskSet)
                    selectorCell = ccObj.getSelectorFromMasks;
                else
                    for i=1:ccObj.numCov
                        %selectorCell{i} = 1:ccObj.covArray{i}.dimension;
                        selectorCell{i} = 1:ccObj.getCov(i).dimension;
                    end
                end
            end
            if(nargin<2)
                repType='standard';
            end            
            
            dimTot = sumDimensions(selectorCell);
            nCov   = numActCov(selectorCell);
            covInd = covIndFromSelector(selectorCell);
            
            dataMat=zeros(length(ccObj.getCov(1).getSigRep.time),dimTot);
%             size(dataMat)
            for i=1:nCov
                if(i==1)
                    currentOffset =0;
                else
                    currentOffset = sumDimensions(selectorCell,covInd(i-1));
                end
                    %covariate.getCovMatrix(covObj,repType, selectorArray,binwidth,minTime,maxTime)
                    data=ccObj.getCov(covInd(i)).getSigRep(repType).dataToMatrix(selectorCell{covInd(i)});%,binwidth,minTime,maxTime);
                    endInd = min(size(dataMat,1),size(data,1));
                    dataMat(1:endInd,currentOffset+(1:length(selectorCell{covInd(i)})))=data(1:endInd,:);
            end
        end
        function structure=dataToStructure(ccObj,selectorCell,binwidth, minTime, maxTime)
            % structure=dataToStructure(ccObj,selectorCell,binwidth, minTime, maxTime)
            % structure representation of the CovColl.
            
            if(nargin<5)
                maxTime  = ccObj.maxTime;
            end
            if(nargin<4)
                minTime = ccObj.minTime;
            end
            if(nargin<3) 
                binwidth = 1/ccObj.getCov(1).getSigRep.sampleRate;
            end
            if(nargin<2)
                if(ccObj.isCovMaskSet)
                   selectorCell = ccObj.getSelectorFromMasks;
                else
                    for i=1:ccObj.numCov
                        %selectorCell{i} = 1:ccObj.covArray{i}.dimension;
                         selectorCell{i} = 1:ccObj.getCov(i).dimension;
                    end
                end
            end
            repType = 'standard';
            dataMatrix =ccObj.dataToMatrix(repType, selectorCell, binwidth,minTime,maxTime);
           
            %Convert to a standard matlab structure
            structure.time=ccObj.getCov(1).time;
            structure.signals.values=dataMatrix;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotHandle = plot(ccObj,handle,repType,selectorCell)
            if(nargin<4)
                if(ccObj.isCovMaskSet)
                    selectorCell = ccObj.getSelectorFromMasks;
                else
                    for i=1:ccObj.numCov
                       %selectorCell{i} = 1:ccObj.covArray{i}.dimension;
                       selectorCell{i} = 1:ccObj.getCov(i).dimension;
                    end
                end
            end                
            if(nargin<3)
                repType='standard';
            end
            if(nargin<2)
                handle = gcf;
            end
            
            %plotHandle = figure(handle);
            nCov=numActCov(selectorCell);
            covInd = covIndFromSelector(selectorCell);
            if(handle == gcf) %given a figure to plot in;
                for i=1:nCov
                    if(nCov==1)
                       %no subplot
                    elseif(nCov==2)
                        subplot(2,1,i)
                    elseif(nCov==3)
                        subplot(3,1,i)
                    elseif(nCov==4)
                        subplot(2,2,i)
                    else 
                        figure;
                    end
                    ch=gca;
                    % h=plot(sObj,selectorArray,plotProps,handle)
                    currentObj = ccObj.getCov(covInd(i));
                    plotHandle=currentObj.plot([],[],ch); %default selectorArray and default plotProps
                end
            elseif(length(handle)==nCov) %got a subplot for each covariate
                for i = 1:length(handle)
                    currentObj = ccObj.getCov(covInd(i)).getSigRep(repType);
                    axes(handle(i));
                    plotHandle=currentObj.plot(selectorCell{i},[],handle(i));
                end
            end
        end
        
       
    end
    
    methods (Access = private)
        function setMasksFromSelector(ccObj,selectorCell)
            if(length(selectorCell)==ccObj.numCov);
                ccObj.covMask=ccObj.getCovMaskFromSelector(selectorCell);
            end
        end        
        function cMask =getCovMaskFromSelector(ccObj,selectorCell)
            cMask = cell(1,ccObj.numCov);
            for i=1:length(cMask)
                cMask{i}=zeros(1,length(ccObj.getCov(i).dataMask));
                if(~isempty(selectorCell{i}))
                    if(length(selectorCell{i})>1  && max(selectorCell{i})==1)
                        cMask{i}(selectorCell{i}==1)=1;
                    else
                        cMask{i}(selectorCell{i})=1;
                    end
                end
            end
        end        
        function selectorArray = getSelectorFromMasks(ccObj,covMask)
            if(nargin<2)
                covMask=ccObj.covMask;
            end
               
            selectorArray=cell(1,ccObj.numCov);
            for i=1:ccObj.numCov
                ind=find(covMask{i}==1);
                if(~isempty(ind))
                   selectorArray{i} = ind;
                else
                    selectorArray{i} = [];
                end
            end
        end
   
        function answer=isaSelectorCell(ccObj,dataSelector)
            if(length(dataSelector)==ccObj.numCov && ~containsChars(dataSelector))
                answer=1;
            else
                answer=0;
            end
        end
        function selectorCell = generateSelectorCell(ccObj, dataSelector)
            %dataSelector must be in the following format
            %dataSelector{1} = {'Position','x','y'};
            %dataSelector{2} = {'Force','fx','fy','fz'};
            selectorCell=cell(1,ccObj.numCov);
            if(isempty(dataSelector))
                for i = 1:length(selectorCell);
                    selectorCell{i} =[]; %zeros(1,ccObj.getCov(i).dimension);
                end
            else
                if(isa(dataSelector{1},'char'))
                    covName=dataSelector{1};
                    covLabels=cell(1,length(dataSelector)-1);
                    for i =1:length(covLabels)
                        covLabels{i}=dataSelector{i+1};
                    end
                    covIndex=ccObj.getCovIndFromName(covName);
                    currCov = ccObj.getCov(covIndex);
                    selectorCell{covIndex}=currCov.getIndicesFromLabels(covLabels);
                elseif(isa(dataSelector{1},'cell'))
                    for i=1:length(dataSelector)
                        [covName, covLabels] = parseDataSelectorArray(dataSelector{i});
                        covIndex=ccObj.getCovIndFromName(covName);
                        currCov = ccObj.getCov(covIndex);
                        if(~isempty(currCov))
                            selectorCell{covIndex}=currCov.getIndicesFromLabels(covLabels);
                        else
                           error(['Covariate ' covName ' not found!']);
                        end
                    end
                elseif(isa(dataSelector{1},'double'))
                    selectorCell=dataSelector;
                else
                    error('dataSelector specified incorrectly!');
                end
            end
            
        end
        function addCovCellToColl(ccObj,cov)
            [~, ncolumns]=size(cov);
                for i=1:ncolumns
                    if(isa(cov{i},'Covariate'))
                        ccObj.addSingleCovToColl(cov{i});
                    else
                        error('CovColl requires a cell array of Covariate class elements');
                    end
                end
        end
        function addSingleCovToColl(ccObj,cov)
            if(~ccObj.isCovPresent(cov))
                ccObj.covArray{ccObj.numCov+1}= cov;
                ccObj.updateTimes(cov);
                ccObj.covDimensions(ccObj.numCov+1) = cov.dimension;
                ccObj.covMask{ccObj.numCov+1} = cov.dataMask;
                ccObj.numCov = ccObj.numCov + 1;
                %ccObj.sampleRate
                %cov.sampleRate
                if(isempty(ccObj.sampleRate)) %this is our first element
                    ccObj.sampleRate = cov.sampleRate;
                    ccObj.originalSampleRate = ccObj.sampleRate;
                elseif(ccObj.sampleRate==cov.sampleRate)
                    %Do nothing - just add
                elseif(ccObj.sampleRate>cov.sampleRate)  %Upsample Covariate
                    cov.setSampleRate(ccObj.sampleRate);
                elseif(ccObj.sampleRate<cov.sampleRate); %Upsample other covariates in collection
                    ccObj.setSampleRate(cov.sampleRate);
                else
                    error('Problem setting the sample rate during adding covariate to collection');
                end
            else
                error('Covariate not added because it is already present in this collection or another covariate has the same name');
            end        
        end
        
        function updateTimes(ccObj,cov)
            
            timeVec=cov.getSigRep.getTime;
            minTime=min(timeVec); maxTime=max(timeVec);
            if(minTime<ccObj.minTime)
                ccObj.setMinTime(minTime);
            end
            if(maxTime>ccObj.maxTime)
                ccObj.setMaxTime(maxTime);
            end
        end
        function ind = getCovIndFromName(ccObj,name)
            ind=[];
            for i=1:ccObj.numCov
                if(strcmp(ccObj.getCov(i).name,name))
                    ind=i;
                    break;
                end
            end
        end
          function removeFromColl(ccObj,identifier)
                covs = ccObj.getCov(identifier);
                ind = zeros(1,length(covs));
                if(length(ind)>1)
                    for i = 1:length(ind)
                        ind(i) = ccObj.getCovIndFromName(covs{i}.name);
                    end
                else
                    ind=ccObj.getCovIndFromName(covs.name);
                end
                
                ccObj.removeFromCollByIndices(ind);
                
        end        
        function removeFromCollByIndices(ccObj,ind)
            remaining = ccObj.generateRemainingIndex(ind);
            covArray = cell(1,length(remaining));
            covMask = cell(1,length(remaining));
            covDimensions = zeros(1,length(remaining));
            for i=1:length(remaining)
                cov = ccObj.getCov(remaining(i));
                covMask{i} = ccObj.covMask{remaining(i)};
                covArray{i} = cov;
                covDimensions(i) = cov.dimension;
            end
            numCov = length(remaining);
            ccObj.covArray = covArray;
            ccObj.covMask  = covMask;
            ccObj.numCov = numCov;
            ccObj.covDimensions = covDimensions;
            minTime=ccObj.findMinTime;
            maxTime=ccObj.findMaxTime;
            ccObj.setMinTime(minTime);
            ccObj.setMaxTime(maxTime);
            if(numCov==0)
                ccObj.sampleRate =[];
                ccObj.originalSampleRate = [];
            end
        
        end
        function remain = generateRemainingIndex(ccObj,ind)
            remain=zeros(1,ccObj.numCov-length(ind));
            count=1;
            for i=1:ccObj.numCov
                if(sum(i==ind)>0) %then this is one of the indices we are removing
                 % do nothing
                else
                    remain(count) = i;
                    count=count+1;
                end
            end 
        end
    end
    methods (Static)
        function ccObj = fromStructure(structure)
            if(isa(structure,'struct'))
                cov = cell(1,structure.numCov);
                for i=1:structure.numCov;
                    cov{i} = Covariate.fromStructure(structure.covArray{i});
                end
                ccObj = CovColl(cov);
    %             covMask = structure.covMask;
    %             ccObj.setMask(covMask); 
                %% Need to fix how mask is set!!!
                ccObj.setMinTime(structure.minTime);
                ccObj.setMaxTime(structure.maxTime);
            elseif(isa(structure,'cell'))
                ccObj= cell(length(structure),1);
               for i=1:length(structure)
                  ccObj{i} = CovColl.fromStructure(structure{i}); 
               end
                
            end
        end
        
    end
    
    
    
end

%Helper functions

function ind = covIndFromSelector(selectorCell)
  ind=zeros(1,numActCov(selectorCell));
  count=1;
  for i=1:length(selectorCell)
      if(~isempty(selectorCell{i}))
          ind(count)=i;
          count=count+1;
      end
  end
end
function n = numActCov(selectorCell)
    n=0;
    for i=1:length(selectorCell)
        if(~isempty(selectorCell{i}))
            n=n+1;
        end
    end
end

function dimTot = sumDimensions(selectorCell,index)
    if(nargin<2)
        index=length(selectorCell);
    end
    
    dimTot=0;
    if(index>0 && index<=length(selectorCell))
        for i=1:index
            dimTot=dimTot+length(selectorCell{i});
        end
    end
end
function [covName, covLabels] = parseDataSelectorArray(entry)
    covName = entry{1};
    covLabels = cell(1,length(entry)-1);
    for i =1:length(covLabels)
        covLabels{i} = entry{i+1};
    end
end
function answer=containsChars(x)
    if(isa(x,'cell'))
        for i=1:length(x)
            if(isa(x{i},'char'))
                answer=1;
                break;
            end
        end
        answer=0;
    elseif(isa(x,'char'))
        answer =1;
    else
        answer=0;
    end
 end

    
    


