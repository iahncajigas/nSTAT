classdef FitResult < handle
% FITRESULT 
% stores results of a fit using the Analysis object
% The results are for a single neuron over a range of configurations
%
% <a href="matlab: methods('FitResult')">methods</a>
% <a href="matlab:web('FitResultExamples.html', '-helpbrowser')">FitResult Examples</a> 
%
% see also <a href="matlab:help('Analysis')">Analysis</a>
%
% Reference page in Help browser
% <a href="matlab:doc('FitResult')">doc FitResult</a>


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
       numResults   %Number of results in this FitResult object
       lambda       %Lambda signal
       numCoeffs    %Number of coefficients for each fitResult
       fitType      %Poisson or Binomial
       
       b            %coefficients for each fit
       
       dev          %deviance for each fit
       AIC          %Akaike's Information Criterion for each fit
       BIC          %Baysian Information Criterion for each fit
       stats        %Relevant statistics for each fit
       configs      % the config collection for the different fits
       configNames  % names of the the differen fits
       neuronNumber % the number of the neuron the data comes from
       neuralSpikeTrain % the spike data
       covLabels    
       uniqueCovLabels
       indicesToUniqueLabels
       numHist        % Number of history terms (used for indexing into the regression coefficients) for each fit
       histObjects    % History object for self firing 
       ensHistObjects % History object to be applied to the neuron's neigbors to compute ensemble effect
       flatMask
       Z    % Rescaled spike times from the Time-Rescaling theorem, exponential rate 1
       U    % Transformed z's -> uniform in [0,1]
       X    % Transformed u's -> gaussian in [-inf, inf]
       Residual %fit residual
%        xAxis
%        KSSorted
%        ks_stat
       invGausStats 
       KSStats  %Kolmogorov Smirnov Statistics
       plotParams
       XvalData % cell array of raw data used for validation
       XvalTime % cell array of time vectors for each element of XvalData
       validation % a FitResult object with the validation data
       minTime  % The minTime from the spikeTrain (not necessarily the analysis)
       maxTime  % The maxTime for the spikeTrain (not necessarily the analysis)
    end
    properties (Constant,Hidden)
        colors={'b','g','r','c','m','y','k'};    
    end


    methods 
     
        function fitObj=FitResult(spikeObj,covLabels,numHist,histObjects,ensHistObj,lambda,b, dev, stats,AIC,BIC,configColl,XvalData,XvalTime,distribution)
            % fitObj=FitResult(spikeObj,covLabels,numHist,histObjects,ensHistObj,lambda,b, dev, stats,AIC,BIC,configColl,XvalData,XvalTime)
            % Stores the results of multiple regressions for a single  neuron into a accessible structure.
            %
            % spikeObj: The spike train for the neuron whose results are
            %           being stored.
            % covLabels: A 2-d cell array, the jth row has all the labels for the covariates used in the jth fit
            % numHist: The number of history terms in each of the N fits.
            % histObjects: The History object for each of the N fits
            % ensHistObj: The History object for used to compute the ensemble history effect.
            % lambda: The conditional intensity function evaluated usin the
            % data. Each dimension of lambda corresponds to the a different
            %       GLM Fit.
            % b: N-component cell array containing the GLM regression
            %    coefficient of each of the fits. The jth component has all
            %    the regression coefficients for the jth trial.
            % dev: vector of Deviances for each the GLM fits.
            % stats: Cell array of the stats parameters for each GLM fit;
            % AIC: vector of Akaike's information criteria for each the GLM fits.
            % BIC: vector of Bayes Information criteria for each the GLM fits.
            % configColl: configCollection object used to generate this
            %             results
            % XvalData: Data to be used for validation.
            % XvalTime: Time vector for the data.
            
            if(nargin< 14)
                XvalTime =[];
            end
            if(nargin<13)
                XvalData =[];
            end
            
            if(isa(spikeObj,'cell'))
                for i=1:length(spikeObj)
                    if(isnumeric(spikeObj{i}.name))
                        nNumber(i) =spikeObj{i}.name;
                    else
                        nNumber(i) = str2double(spikeObj{i}.name(~isletter(spikeObj{i}.name)));
                    end
                    minTime(i)=spikeObj{i}.minTime;
                    maxTime(i)=spikeObj{i}.maxTime;
                end
                nNumber = unique(nNumber);
                minTime = unique(minTime);
                maxTime = unique(maxTime);
                if(length(nNumber)>1)
                    error('Can only have a FitResults with spike trains from a single neuron');
                end
                if(length(minTime)>1 || length(maxTime)>1)
                    error('Spike Trains are of different lengths');
                end
                
            elseif(isa(spikeObj,'nspikeTrain'))
                if(isnumeric(spikeObj.name))
                    nNumber =spikeObj.name;
                else
                    nNumber = str2double(spikeObj.name(~isletter(spikeObj.name)));
                end
                minTime=spikeObj.minTime;
                maxTime=spikeObj.maxTime;
                    
            end
                
            fitObj.neuronNumber = nNumber; %str2num(spikeObj.name);
            fitObj.neuralSpikeTrain = spikeObj;
            fitObj.minTime = minTime;
            fitObj.maxTime = maxTime;

            fitObj.numResults = 0;
            fitObj.configs = configColl;
            fitObj.configNames = configColl.getConfigNames;
            fitObj.covLabels=covLabels;
            fitObj.uniqueCovLabels= getUniqueLabels(covLabels);
           
            fitObj.mapCovLabelsToUniqueLabels;
            fitObj.numHist=numHist;
            fitObj.histObjects = histObjects;
            fitObj.ensHistObjects = ensHistObj;
            fitObj.addParamsToFit(fitObj.neuronNumber,lambda,b, dev, stats,AIC,BIC,configColl);
            fitObj.Z        =[]; %rescaled spikes times - exponentially dist.
            fitObj.U        =[]; %rescaled spike times - uniformly dist.
            fitObj.X        =[]; %rescaled spike times - gaussian dist.
            fitObj.Residual =[]; %fit residual for PP
            fitObj.KSStats.xAxis    =[];
            fitObj.KSStats.KSSorted =[];
            fitObj.KSStats.ks_stat  =[];
            fitObj.invGausStats.rhoSig=[];
            fitObj.invGausStats.confBoundSig=[];
            fitObj.plotParams = [];
            fitObj.XvalData = XvalData;
            fitObj.XvalTime = XvalTime;     
            
            fitObj.fitType = distribution;     
             
        end       
        function fitObj = setNeuronName(fitObj,name)
           fitObj.neuronNumber = name; 
        end
        function mFitRes = mergeResults(fitObj,newFitObj)
            % mFitRes = mergeResults(fitObj,newFitObj)
            % mFitRes contains the results from fitObj followed by the
            % results of newFitObj in a single new FitResult object
            %
            % newFitObj can be of class 'FitResult' or a cell array of
            % 'FitResult' objects. In the latter case, the results are
            % apppended in the order that they appear in each cell array.
            if(isa(newFitObj,'FitResult'))
%                 newFitObj.neuronNumber
                if(fitObj.neuronNumber ==newFitObj.neuronNumber)
                    spikeObj = fitObj.neuralSpikeTrain;
                    covLabels = fitObj.covLabels(1:fitObj.numResults);
                    covLabels((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.covLabels(1:newFitObj.numResults);
                    numHist = fitObj.numHist(1:fitObj.numResults);
                    numHist((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.numHist(1:newFitObj.numResults);
                    histObjects=fitObj.histObjects(1:fitObj.numResults);
                    histObjects((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.histObjects(1:newFitObj.numResults);
                    ensHistObjects=fitObj.ensHistObjects(1:fitObj.numResults);
                    ensHistObjects((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.ensHistObjects(1:newFitObj.numResults);
                    b=fitObj.b(1:fitObj.numResults);
                    b((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.b(1:newFitObj.numResults);
                    dev = [fitObj.dev newFitObj.dev];
                    AIC = [fitObj.AIC newFitObj.AIC];
                    BIC = [fitObj.BIC newFitObj.BIC];
                    stats=fitObj.stats(1:fitObj.numResults);
                    stats((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.stats(1:newFitObj.numResults);
                    lambda = fitObj.lambda.merge(newFitObj.lambda);
                    
                    for i=1:fitObj.numResults
                        config{i}=fitObj.configs.getConfig(i);
                    end
                    offset=fitObj.numResults;
                    for i=1:newFitObj.numResults
                        config{i+offset}=newFitObj.configs.getConfig(i);
                    end
                    configColl= ConfigColl(config);
                    
                    XvalData = [fitObj.XvalData newFitObj.XvalData];
                    XvalTime = [fitObj.XvalTime newFitObj.XvalTime];
                    distribution=fitObj.fitType(1:fitObj.numResults);
                    distribution((fitObj.numResults+1):(fitObj.numResults+newFitObj.numResults)) = newFitObj.fitType(1:newFitObj.numResults);
                    tempZ = zeros(length(fitObj.Z),size(newFitObj.Z,2));
                    tempU = zeros(length(fitObj.U),size(newFitObj.U,2));
                    tempZ(1:length(newFitObj.Z),:) = newFitObj.Z;
                    tempU(1:length(newFitObj.U),:) = newFitObj.U;
                    Z=[fitObj.Z tempZ];
                    U=[fitObj.U tempU];
                    [X,rhoSig,confBoundSig] = Analysis.computeInvGausTrans(Z);

                    M=fitObj.Residual.merge(newFitObj.Residual);
                   
                    origLength = size(fitObj.KSStats.xAxis,1);
                    currLength = size(newFitObj.KSStats.xAxis,1);
                    
                    if(currLength~=origLength)
                        %we use this because some times the time scales
                        %dont match up.  In particular when spikeTrain is
                        %segmented by steps or windows and the window sizes
                        %are normalized to 1.
                      newX = fitObj.KSStats.xAxis;
                      oldX = newFitObj.KSStats.xAxis;
                      oldY = newFitObj.KSStats.KSSorted;
                      y = interp1(oldX,oldY,newX(:,1),'spline','extrap');
                      xAxis    = [fitObj.KSStats.xAxis newX(:,1)];
                      KSSorted = [fitObj.KSStats.KSSorted y];
                    else
                       
                        xAxis    = [fitObj.KSStats.xAxis newFitObj.KSStats.xAxis];
                        KSSorted = [fitObj.KSStats.KSSorted newFitObj.KSStats.KSSorted];
                        
                    end
                    ks_stat  = [fitObj.KSStats.ks_stat newFitObj.KSStats.ks_stat];
                    mFitRes=FitResult(spikeObj,covLabels,numHist,histObjects,ensHistObjects,lambda,b, dev, stats,AIC,BIC,configColl,XvalData,XvalTime,distribution);
                    mFitRes.setKSStats(Z,U, xAxis, KSSorted, ks_stat);
                    mFitRes.setInvGausStats(X,rhoSig,confBoundSig);
                    mFitRes.setFitResidual(M);
                    
                    
                elseif(isa(newFitObj,'cell'))
                    if(isa(newFitObj{1},'FitResult'))
                        for i=1:length(newFitObj)
                            if(i==1)
                                mFitRes = fitObj.mergeResults(newFitObj{i});
                            else
                                mFitRes = mFitRes.mergeResults(newFitObj{i});
                            end
                        end
                    end
                end
            end
            
        end
        function addParamsToFit(fitObj,neuronNum,lambda,b, dev, stats,AIC,BIC,configColl)
            % addParamsToFit(fitObj,neuronNum,lambda,b, dev, stats,AIC,BIC,configColl)
            % Add the specified parameters to the current FitResult object
            % only if the neuronNum matches the neuronNum of this object
            if(fitObj.neuronNumber==neuronNum)
              if(isa(lambda,'cell'))
                  newLambda=lambda{1};
                  for i=2:length(lambda)
                      newLambda = newLambda.merge(lambda{i});
                  end
              elseif(isa(lambda,'Covariate')||isa(lambda,'SignalObj'))
                  newLambda = lambda;
              end

              numNewResults = newLambda.dimension;%number of new elements
              if(nargin<8)
                  configColl = cell(1,numNewResults);
              end

              if(numNewResults==1)
                      fitObj.b{fitObj.numResults+1}    = b{1};
                      
                      fitObj.dev(fitObj.numResults+1)  = dev;
                      fitObj.stats{fitObj.numResults+1}= stats{1};
                      if(nargin<7)
                          fitObj.AIC(fitObj.numResults+1)  = 2*length(b)+dev;
                          fitObj.BIC(fitObj.numResults+1)  = length(b)*log(length(newLambda.time))+dev;
                      else
                          fitObj.AIC(fitObj.numResults+1)  = AIC;
                          fitObj.BIC(fitObj.numResults+1)  = BIC;
                      end
                          
                      fitObj.numCoeffs(fitObj.numResults+1) = length(b);
              else
                  for i=1:numNewResults
                      fitObj.b{fitObj.numResults+i}    = b{i};
                      
                      fitObj.dev(fitObj.numResults+i)  = dev(i);
                      fitObj.stats{fitObj.numResults+i}= stats{i};
                      if(nargin<7)
                          fitObj.AIC(fitObj.numResults+i)  = 2*length(b{i})+dev(i);
                          fitObj.BIC(fitObj.numResults+i)  = length(b{i})*log(length(newLambda.time))+dev(i);
                      else
                          fitObj.AIC(fitObj.numResults+i)  = AIC(i);
                          fitObj.BIC(fitObj.numResults+i)  = BIC(i);
                      end
                      fitObj.numCoeffs(fitObj.numResults+i) = length(b{i});
                  end
              end
              if(fitObj.numResults ==0)
                  fitObj.lambda = newLambda;
              else
                fitObj.lambda = fitObj.lambda.merge(newLambda); %new lambda
              end
              
              fitObj.numResults = fitObj.numResults+numNewResults;
              dataLabels = cell(1,fitObj.numResults);
              for i=1:fitObj.numResults
                 dataLabels{i} = strcat('\lambda_{',num2str(i),'}');
              end
              fitObj.lambda.setDataLabels(dataLabels);
              
              fitObj.configs.addConfig(configColl);
              fitObj.configNames = fitObj.configs.getConfigNames;
          else
              error('Neuron number does not match');
          end
        end
        function lambda = computeValLambda(fitObj)
            % lambda = computeValLambda(fitObj)
            % Returns a Covariate object lambda. This is the Conditional
            % intensity function evaluated using the validation data
            lambdaData = zeros(length(fitObj.XvalTime{1}),fitObj.numResults);
            for i=1:fitObj.numResults
                lambdaData(:,i) = fitObj.evalLambda(i,fitObj.XvalData{i});
            end
            lambda=Covariate(fitObj.XvalTime{1},lambdaData,...
                  '\lambda(t)',fitObj.lambda.xlabelval,...
                  fitObj.lambda.xunits,'Hz',fitObj.lambda.dataLabels);
        end
        
        function mapCovLabelsToUniqueLabels(fitObj)
            % mapCovLabelsToUniqueLabels(fitObj)
            % Used internally by the FitResult class generate a matrix that
            % maps how covariate labels of the fit object map to unique
            % covariate labels. For example, multiple fits that have a
            % constant baseline term will be assumed to refer to the same
            % "baseline" term and not two separate ones
            flatMask = zeros(length(fitObj.uniqueCovLabels),length(fitObj.covLabels));
            for j=1:length(fitObj.covLabels)
                currLabels = fitObj.covLabels{j};
                index=zeros(1,length(currLabels));
                for i=1:length(currLabels)
                    index(i)=strmatch(currLabels{i}, fitObj.uniqueCovLabels, 'exact');
                end
                
                fitObj.indicesToUniqueLabels{j} = index;
                flatMask(index,j) = 1;
            end
            fitObj.flatMask = flatMask;
        end        
        function p=getPlotParams(fitObj)
            % p=getPlotParams(fitObj)
            if(isempty(fitObj.plotParams))
                fitObj.computePlotParams;
            end
                p=fitObj.plotParams;
        end        
        function plotValidation(fitObj) 
            % plotValidation(fitObj) 
            % calls plotResults on the validation FitResult object if
            % validation data is present. Note that the GLM coefficients
            % are not recomputed and therefore the same as those obtained
            % from the training data.
            if(~isempty(fitObj.validation))
                fitObj.validation.plotResults;
            else
                display('Validation Data not available to plot');
            end
        end
        function answer = isValDataPresent(fitObj)
            % answer = isValDataPresent(fitObj)
            % returns 1 if validation data is present. This method is used
            % to determine if validation data is available to compute the
            % validation results.
            answer = 0;
            if(~isempty(fitObj.XvalTime) && ~isempty(fitObj.XvalData))
                for i=1:length(fitObj.XvalTime)
                    currTime = fitObj.XvalTime{i};
                    if(~isempty(currTime))
                        if(currTime(end)-currTime(1)>0)
                            answer =1;
                            break;
                        end
                    end
                end
                
            end
            
        end        
        function lambdaData = evalLambda(fitObj,lambdaIndex,newData)
            % lambdaData = evalLambda(fitObj,lambdaIndex,newData)
            % lambdaIndex: the index of the corresponding lambda to be
            %              evaluated with the new data.
            % newData:     matrix of covariates in same order as fits without
            %              constant term in first column
%             if(isa(newData,'double'))
%                 [~,columns] = size(newData);
%                 tempData = cell(1,columns);
%                 for i=1:columns
%                    tempData{i} = newData(:,i); 
%                 end
%                 newData = tempData;
%             end
            
            if(lambdaIndex>0 && lambdaIndex <= fitObj.numResults)
                b=fitObj.b{lambdaIndex}; %coefficient matrix
                if(isempty(newData))
                    [rows,~] = size(newData);
                    baseline=ones(rows,1);
                    lambdaData =  exp(b(1)*baseline);
                else
                    if(isa(newData,'double')) %matrix, 1 column per coefficient
                        baseline=ones(length(newData),1);
                        [~,columns] = size(newData);
                        
                        if(length(b)>=1)
                            lambdaData = exp(newData*b(1:end));
                            if(strcmp(fitObj.fitType{lambdaIndex},'poisson'))
%                                 lambdaData = exp(newData*b(1:end));
%                                 lambdaData = exp(b(1) + newData*b(2:end));
                            else 
%                                 lambdaData = exp(b(1) + newData*b(2:end));
                                lambdaData = lambdaData./(1+lambdaData);
                            end
%                         else
%                             if(strcmp(fitObj.fitType{lambdaIndex},'poisson'))
%                                 lambdaData = exp(b(1)*baseline);
%                                 
%                             else
%                                 lambdaData = exp(b(1)*baseline);
%                                 lambdaData = lambdaData./(1+lambdaData);
%                             end
                        end
                        lambdaData = lambdaData*fitObj.neuralSpikeTrain.sampleRate;
                    elseif(isa(newData,'cell')) % a cell array, each element is matrix of values for each coeff
%                         baseline=ones(size(newData{1})); %design matrix
                            runSum=0;
                        for i=1:(length(newData)) %-fitObj.numHist(lambdaIndex))
%                             if(i==1)
%                                 runSum = b(1)*baseline;
%                             else
                                if(i<=length(b))
                                    runSum = runSum+b(i)*newData{i};
                                end
%                             end
                        end
                        if(strcmp(fitObj.fitType{lambdaIndex},'poisson'))
                            lambdaData = exp(runSum);
                            lambdaData = lambdaData*fitObj.neuralSpikeTrain.sampleRate;
                        else
                            lambdaData = exp(runSum);
                            lambdaData = lambdaData./(1+lambdaData);
                            lambdaData = lambdaData*fitObj.neuralSpikeTrain.sampleRate;
                        end
                    else
                        error('New data must be cell or a matrix');
                    end
                    
                end
            else
                error('Index into fit params is incorrect');
            end
            
        end        
%         function handle = plotHist(fitObj,fitNum)
%            % handle = plotHist(fitObj,fitNum)
%            % plots the history terms used in this FitResult object
%            % if fitNum is not specified then fitNum=1:numResults
%            if(nargin<2 || isempty(fitNum))
%                fitNum = 1:fitObj.numResults;
%            end
%            
%            for j=fitNum
%                if(j>0 && j <= fitObj.numResults)
%                     b=fitObj.b{j}; %coefficient matrix
%                     startHistIndex = length(b)-fitObj.numHist(j)+1;
%                     if(startHistIndex<length(b))
%                         bHist = b(startHistIndex:end);
%                         if(~isempty(fitObj.histObjects{j}))
%                             windowTimes = fitObj.histObjects{j}.windowTimes;
%                             t=linspace(windowTimes(1),windowTimes(end),100)';
%                             histEffect = zeros(length(t),1);
%                             for i=1:length(windowTimes)-1
%                                 index = and(t>=windowTimes(i),t<=windowTimes(i+1));
%                                 histEffect(index)=exp(-bHist(i))-1; %To offset zero coeffs
%                             end
%                         end
%                     else
%                         t=[0; 0.00001];
%                         histEffect =[0;0];
%                     end
%                     if(j==fitNum(1))
%                         hSig = SignalObj(t,histEffect,'History','time','s','',fitObj.lambda.dataLabels{j});
%                     else
%                         hSig = hSig.merge(SignalObj(t,histEffect,'History Effect','time','s','',fitObj.lambda.dataLabels{j}));
%                     end
%                end
%            end
%            N=floor(length(hSig.time)./70); B=ones(1,N)/N; A=1;
%            handle=hSig.filtfilt(B,A).plot;    
%         end
        function computePlotParams(fitObj,fitNum)
             if(nargin<2)
               fitNum = 1:fitObj.numResults;
             end
           index=find(sum(fitObj.flatMask,2)>0);%1:length(fitObj.flatMask(:,1));
           %Only use the labels that appear in at least one fit
           %Otherwise that parameter was not present for any of the
           %regressions and just takes up plot real-estate
           
           sigIndex=zeros(length(index),length(fitNum));
           bAct = nan(length(index),length(fitNum));
           seAct= nan(length(index),length(fitNum));
           
           for i=fitNum
               %this indexing is to avoid extremely large se's from
               %affecting plots
               criteria = find(fitObj.stats{i}.se'<100);
               %indicesForFit = find(fitObj.flatMask(index,i)==1);
               indicesForFit = fitObj.indicesToUniqueLabels{i};
               bVals = fitObj.b{i}(criteria);
               bAct(indicesForFit(criteria),i) = bVals; %sorted according to uniqueLabels
               seVals = fitObj.stats{i}.se(criteria)';
               seAct(indicesForFit(criteria),i)= seVals; %sorted according to uniqueLabels;
               temp = sign([bAct(:,i)-seAct(:,i) bAct(:,i)+seAct(:,i)]);
               productOfSigns = temp(:,1).*temp(:,2); %should be positive
               sIndex=and(productOfSigns>0,seAct(:,i)~=0);
               sigIndex(:,i)=sIndex;
           end
           fitObj.plotParams.bAct = bAct;
           fitObj.plotParams.seAct= seAct;
           fitObj.plotParams.sigIndex = sigIndex;
           fitObj.plotParams.xLabels  = cell(length(index),1);           
           fitObj.plotParams.xLabels = fitObj.uniqueCovLabels;
               
%            for i=1:(length(index))
%                if(i==1)
%                    fitObj.plotParams.xLabels{i} = 'baseline';
%                    %text(i, 0,'baseline','interpreter','latex');
%                else
%                    fitObj.plotParams.xLabels{i} = fitObj.covLabels{index(i)-1};
%                    %text(i, 0,fitObj.covLabels{index(i)-1},'interpreter','latex');
%                end
%            end
           tempVal =sum(fitObj.flatMask,2);
           fitObj.plotParams.numResultsCoeffPresent =tempVal(index); 
        end
        
        function [coeffIndex, epochId,numEpochs] = getCoeffIndex(fitObj,fitNum,sortByEpoch)
          if(nargin<3 || isempty(sortByEpoch))
             sortByEpoch=0; 
          end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:fitObj.numResults;
          end
          if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
          end   
          [histIndex, epochId] = fitObj.getHistIndex(fitNum,sortByEpoch);
          allIndex = 1:length(fitObj.uniqueCovLabels);
          
          nonHistIndex = setdiff(allIndex,histIndex);
%           nonNANIndex = find(sum(~isnan(fitObj.plotParams.bAct(:,fitNum)),2)>=1);
          nonNANIndex=  allIndex;
          actCoeffIndex = nonHistIndex(ismember(nonHistIndex, nonNANIndex));
          allCoeffTerms = fitObj.uniqueCovLabels(actCoeffIndex);
%           coeffName = cell(size(allCoeffTerms));
          epochStartInd=regexp(allCoeffTerms,'_\{\d*\}','start'); 
          epochEndInd=regexp(allCoeffTerms,'_\{\d*\}','end');
         
          allCoeffIndex = [];
          nonEpochIndex=[];
%           nonEmptyCoeffNameInd = [];
          epochsExist =0;
          for i=1:length(allCoeffTerms)
              if(~isempty(allCoeffTerms{i}))
                allCoeffIndex = [allCoeffIndex i];
                
                 if(~isempty(epochStartInd{i}))
%                      nonEmptyCoeffNameInd = [nonEmptyCoeffNameInd i];
                     epochsExist=1;
                     actStart = epochStartInd{i}+2;
                     actEnd   = epochEndInd{i}-1;
                     numEpoch(i) = str2num(allCoeffTerms{i}(actStart:actEnd));
%                      coeffName{i} =  allCoeffTerms{i}(1:actStart-3);
                     
                     
                 else
                    nonEpochIndex = [nonEpochIndex i];
                    numEpoch(i) = 0; % make terms that only appear once part of epoch 0.
                    
                 end
              end
              
          end
        
%           coeffName = coeffName(nonEmptyCoeffNameInd);
          if(epochsExist && ~sortByEpoch)
              totalEpochs = unique(numEpoch);
              coeffIndex = nonEpochIndex;
              if(nargout>1)
                  epochId=zeros(size(nonEpochIndex));
              end
              for i=1:length(totalEpochs)
                  if(totalEpochs(i)~=0)
                    coeffIndex = [coeffIndex, find(numEpoch==totalEpochs(i))];
                  
                     if(nargout>1)
                        epochId = [epochId, totalEpochs(i)*ones(size(find(numEpoch==totalEpochs(i))))];
                     end
                  end
              end
              coeffIndex = actCoeffIndex(coeffIndex);
          elseif(epochsExist && sortByEpoch)
              coeffIndex = actCoeffIndex(allCoeffIndex);
              if(nargout>1)
                epochId = numEpoch;
              end
          else
              coeffIndex = actCoeffIndex(allCoeffIndex);
              if(nargout>1)
                epochId = zeros(size(allCoeffIndex)); %no epochs exist so just create same index for all;
              end
          end
          

%           nonNANIndex = find(sum(~isnan(fitObj.plotParams.bAct(:,fitNum)),2)>=1);
          nonNANIndex = allIndex;
          coeffIndex = coeffIndex(ismember(coeffIndex, nonNANIndex));
            
          if(nargout>2)
              numEpochs = length(unique(epochId));
          end
          
        end
        
        function h=plotCoeffsWithoutHistory(fitObj,fitNum,sortByEpoch,plotSignificance)
           if(nargin<4 || isempty(plotSignificance))
               plotSignificance=1;
           end
           if(nargin<3 || isempty(sortByEpoch))
              sortByEpoch = 0; 
           end
           if(nargin<2 || isempty(fitNum))
               fitNum = 1:fitObj.numResults;
           end
           if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
           end    
               
            
          coeffIndex = fitObj.getCoeffIndex(fitNum,sortByEpoch);
          h=fitObj.plotCoeffs([],fitNum,[],plotSignificance,coeffIndex);
          
          
            
        end
        
        function [histIndex, epochId,numEpochs] = getHistIndex(fitObj,fitNum,sortByEpoch)
            %if sortByEpoch==1 then we group all regression terms with the
            %same name one next to each other by epoch (time interval).
            %Otherwise, we show all epoch one terms, followed by all epoch
            %2 terms, etc.
           if(nargin<3 || isempty(sortByEpoch))
            sortByEpoch = 0;
           end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:fitObj.numResults;
          end
          if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
          end  
          
          
          allHistTerms = regexp(fitObj.uniqueCovLabels,'^[\w*');
          epochStartInd=regexp(fitObj.uniqueCovLabels,'\]_\{\d*\}','start'); 
          epochEndInd=regexp(fitObj.uniqueCovLabels,'\]_\{\d*\}','end');
          allHistIndex = [];
          epochsExist =0;
          for i=1:length(allHistTerms)
              if(~isempty(allHistTerms{i}))
                allHistIndex = [allHistIndex i];
                 if(~isempty(epochStartInd{i}))
                     epochsExist=1;
                     actStart = epochStartInd{i}+3;
                     actEnd   = epochEndInd{i}-1;
                     numEpoch(i) = str2num(fitObj.uniqueCovLabels{i}(actStart:actEnd));
                 end
              end
              
          end
          
          if(epochsExist && ~sortByEpoch)
              totalEpochs = unique(numEpoch);
              histIndex = [];
              if(nargout>1)
                  epochId=[];
              end
              for i=1:length(totalEpochs)
                 histIndex = [histIndex, find(numEpoch==totalEpochs(i))]; 
                 if(nargout>1)
                    epochId = [epochId, totalEpochs(i)*ones(size(find(numEpoch==totalEpochs(i))))];
                 end
              end
          elseif(epochsExist && sortByEpoch)
              histIndex = allHistIndex;
              if(nargout>1)
                epochId = numEpoch;
              end
          else
              histIndex = allHistIndex;
              if(nargout>1)
                epochId = zeros(size(allHistIndex)); %no epochs exist so just create same index for all;
              end
          end
             
                   
%           nonNANIndex = find(sum(~isnan(fitObj.plotParams.bAct(:,fitNum)),2)>=1);
%           histIndex = histIndex(ismember(histIndex, nonNANIndex));

          if(nargout>2)
              numEpochs = length(unique(epochId));
          end
            
        end
        function [coeffMat, labels, SEMat] = getCoeffs(fitObj, fitNum)
             if(nargin<2 || isempty(fitNum))
                fitNum =1:fitObj.numResults;
            end
            sortByEpoch = 0; % Make sure we have different time series if the history is divided into epochs;
            [coeffIndex, epochId, numEpochs] = fitObj.getCoeffIndex(fitNum,sortByEpoch);
            epochNums = unique(epochId);
            
            
             coeffStrings = fitObj.uniqueCovLabels(coeffIndex);
            baseStringEndIndex =regexp(coeffStrings,'_\{\d*\}','start');
            
            for i=1:length(baseStringEndIndex)
                if(~isempty(baseStringEndIndex{i}))
                    baseStrings{i} = coeffStrings{i}(1:baseStringEndIndex{i}-1);
                else
                    baseStrings{i} = coeffStrings{i};
                end
            end
            uniqueCoeffs = unique(baseStrings);
            
            for i=1:length(uniqueCoeffs)
               coeffStrIndex{i} = coeffIndex(strcmp(baseStrings,uniqueCoeffs{i})); 
               if(min(epochId)==0)
                epochIndices{i} = epochId(strcmp(baseStrings,uniqueCoeffs{i}))+1;
               else 
                epochIndices{i} = epochId(strcmp(baseStrings,uniqueCoeffs{i}));
               end
            end
            
%             
%             for i=1:numEpochs
%                 epochIndices{i} = find(epochId==epochNums(i));
%                 epochLength(i) = length(epochIndices{i});
%             end
            
          
            
            coeffIndMat= nan(length(uniqueCoeffs),numEpochs);
            
            labels = cell(size(coeffIndMat));
            for i=1:length(uniqueCoeffs)
               coeffIndMat(i,epochIndices{i}) = coeffStrIndex{i};
               labels(i,epochIndices{i}) = fitObj.uniqueCovLabels(coeffStrIndex{i});
            end
            
            
%             for i=1:numEpochs
%                coeffIndMat(1:epochLength(i),i) = coeffIndex(epochIndices{i});
%                labels(1:epochLength(i),i) = fitObj.uniqueCovLabels(coeffIndMat(1:epochLength(i),i));
%             end
            
            
            
            if(length(fitNum)>1)
                coeffMat = nan(size(coeffIndMat,1),size(coeffIndMat,2), length(fitNum));
                SEMat    = nan(size(coeffIndMat,1),size(coeffIndMat,2), length(fitNum));
                for i=1:length(fitNum)
                    for j=1:length(uniqueCoeffs)
                        bTemp=fitObj.plotParams.bAct(coeffStrIndex{j},i);    
                        seTemp=fitObj.plotParams.seAct(coeffStrIndex{j},i);    
                        coeffMat(j,epochIndices{j},i) = bTemp';
                        SEMat(j,epochIndices{j},i) = seTemp';
                    end
                end
            else
                coeffMat = nan(size(coeffIndMat,1),size(coeffIndMat,2));
                SEMat    = nan(size(coeffIndMat,1),size(coeffIndMat,2));
                for j=1:length(uniqueCoeffs)
                    bTemp=fitObj.plotParams.bAct(coeffStrIndex{j},fitNum); 
                    seTemp = fitObj.plotParams.seAct(coeffStrIndex{j},fitNum); 
                    coeffMat(j,epochIndices{j}) = bTemp';
                    SEMat(j,epochIndices{j}) = seTemp';
                end

            end
            
        end
        
        function [histMat, labels, SEMat] = getHistCoeffs(fitObj,fitNum)
              if(nargin<2 || isempty(fitNum))
                fitNum =1:fitObj.numResults;
            end
            sortByEpoch = 0; % Make sure we have different time series if the history is divided into epochs;
            [histIndex, epochId, numEpochs] = fitObj.getHistIndex(fitNum,sortByEpoch);
            epochNums = unique(epochId);
            
            
            histcoeffStrings = fitObj.uniqueCovLabels(histIndex);
            baseStringEndIndex =regexp(histcoeffStrings,'_\{\d*\}','start');
            baseStrings = cell(length(baseStringEndIndex),1);
            for i=1:length(baseStringEndIndex)
                if(~isempty(baseStringEndIndex{i}))
                    baseStrings{i} = histcoeffStrings{i}(1:baseStringEndIndex{i}-1);
                else
                    baseStrings{i} = histcoeffStrings{i};
                end
            end
            uniqueCoeffs = unique(baseStrings);
            
            for i=1:length(uniqueCoeffs)
               histcoeffStrIndex{i} = histIndex(strcmp(baseStrings,uniqueCoeffs{i})); 
               if(min(epochId)==0)
                epochIndices{i} = epochId(strcmp(baseStrings,uniqueCoeffs{i}))+1;
               else 
                epochIndices{i} = epochId(strcmp(baseStrings,uniqueCoeffs{i}));
               end
            end
            
%             
%             for i=1:numEpochs
%                 epochIndices{i} = find(epochId==epochNums(i));
%                 epochLength(i) = length(epochIndices{i});
%             end
            
          
            
            histcoeffIndMat= nan(length(uniqueCoeffs),numEpochs);
            labels = cell(size(histcoeffIndMat));
%             SEMat = nan(length(uniqueCoeffs),numEpochs);
            for i=1:length(uniqueCoeffs)
               histcoeffIndMat(i,epochIndices{i}) = histcoeffStrIndex{i};
               labels(i,epochIndices{i}) = fitObj.uniqueCovLabels(histcoeffStrIndex{i});
%                SEMat(i,epochIndices{i}) = fitObj.se;
            end
            
            
%             for i=1:numEpochs
%                coeffIndMat(1:epochLength(i),i) = coeffIndex(epochIndices{i});
%                labels(1:epochLength(i),i) = fitObj.uniqueCovLabels(coeffIndMat(1:epochLength(i),i));
%             end
            
            
            
            if(length(fitNum)>1)
                histMat = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2), length(fitNum));
                SEMat   = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2), length(fitNum));
                for i=fitNum
                    for j=1:length(uniqueCoeffs)
                        bTemp=fitObj.plotParams.bAct(histcoeffStrIndex{j},i);    
                        seTemp = fitObj.plotParams.seAct(histcoeffStrIndex{j},i);    
                        histMat(j,epochIndices{j},i) = bTemp';
                        SEMat(j,epochIndices{j},i) = seTemp';
                    end
                end
            else
                histMat = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2));
                SEMat   = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2));
                for j=1:length(uniqueCoeffs)
                    bTemp=fitObj.plotParams.bAct(histcoeffStrIndex{j},fitNum);       
                    seTemp=fitObj.plotParams.seAct(histcoeffStrIndex{j},fitNum);
                    histMat(j,epochIndices{j}) = bTemp';
                    SEMat(j,epochIndices{j}) = seTemp';
                end

            end
            
            
%             
%             if(nargin<2 || isempty(fitNum))
%                 fitNum =1:fitObj.numResults;
%             end
%             sortByEpoch = 0; % Make sure we have different time series if the history is divided into epochs;
%             [histIndex, epochId, numEpochs] = fitObj.getHistIndex(fitNum,sortByEpoch);
%             epochNums = unique(epochId);
%             for i=1:numEpochs
%                 epochIndices{i} = find(epochId==epochNums(i));
%                 epochLength(i) = length(epochIndices{i});
%             end
%             
%             histIndMat= nan(max(epochLength),numEpochs);
%             labels = cell(size(histIndMat));
%             for i=1:numEpochs
%                histIndMat(1:epochLength(i),i) = histIndex(epochIndices{i});
%                labels(1:epochLength(i),i) = fitObj.uniqueCovLabels(histIndMat(1:epochLength(i),i));
%             end
%             
%             
%             
%             if(length(fitNum)>1)
%                 histMat = nan(size(histIndMat,1),size(histIndMat,2), length(fitNum));
%                 for i=1:length(fitNum)
%                     for j=1:numEpochs
%                         bTemp=fitObj.plotParams.bAct(epochIndices{j},i);                
%                         histMat(1:epochLength(j),j,i) = bTemp;
%                     end
%                 end
%             else
%                 histMat = nan(size(histIndMat,1),size(histIndMat,2));
%                 
%                 for j=1:numEpochs
%                     bTemp=fitObj.plotParams.bAct(epochIndices{j},fitNum);                
%                     histMat(1:epochLength(j),j) = bTemp;
%                 end
% 
%             end
            
        end
        function h=plotHistCoeffs(fitObj,fitNum,sortByEpoch,plotSignificance)
          if(nargin<4 || isempty(plotSignificance))
            plotSignificance=1;
          end
          if(nargin<3 || isempty(sortByEpoch))
             sortByEpoch=0; 
          end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:fitObj.numResults;
          end
          if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
          end  
          histIndex = fitObj.getHistIndex(fitNum,sortByEpoch);
          h=fitObj.plotCoeffs([],fitNum,[],plotSignificance,histIndex);
        end
            
            
        function h=plotCoeffs(fitObj,handle,fitNum,plotProps,plotSignificance,subIndex)
           % h=plotCoeffs(fitObj,handle,fitNum,plotProps,plotSignificance)
           % plots the GLM coefficients for each fit along with the
           % confidence intervals. 
           % fitNum: number of the fit to plot. If not specified, all are
           %         plotted.
           % plotProps: properties to use for the making the plot
           % plotSignificance: If 1 then an asterix (*) is place above
           %                   parameters that are statistically different
           %                   from zero with alpha=5%.
           
           
           if(nargin<5 || isempty(plotSignificance))
               plotSignificance = 1;
           end
           
           if(nargin<4 || isempty(plotProps))
               plotProps = [];
           end
           
           if(nargin<3 || isempty(fitNum))
               fitNum = 1:fitObj.numResults;
           end
           
           if(nargin<2 || isempty(handle))
               handle=gca;
           end
           
           if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
           end
           
           if(nargin<6 || isempty(subIndex))
              subIndex = [fitObj.getHistIndex, fitObj.getCoeffIndex];
           end
           
           bAct = fitObj.getPlotParams.bAct(subIndex,fitNum);
           seAct= fitObj.getPlotParams.seAct(subIndex,fitNum);
           sigIndex=fitObj.getPlotParams.sigIndex(subIndex,fitNum);
           
           if(~isempty(plotProps))
               for i=1:length(fitNum)
                    h(i)=errorbar(handle,1:length(subIndex),bAct(:,i),seAct(:,i),plotProps{i}); hold on;
                    set(h(i), 'LineStyle', 'none', 'Marker', '.');%,...
%                         'Linewidth',1,'Marker','o','MarkerSize',6);
                    currColor = get(h(i),'Color');
                    set(h(i),'MarkerEdgeColor',currColor,'MarkerFaceColor',currColor); 
%                     hE= get(h(i),'Children');
%                     errorbarXData = get(hE(2),'XData');
%                     errorbarXData(4:9:end) = errorbarXData(1:9:end) - 0.2;
%                     errorbarXData(7:9:end) = errorbarXData(1:9:end) - 0.2;
%                     errorbarXData(5:9:end) = errorbarXData(1:9:end) + 0.2;
%                     errorbarXData(8:9:end) = errorbarXData(1:9:end) + 0.2;
%                     set(hE(2), 'XData', errorbarXData);
               end
           else
               Xaxis=repmat(1:length(bAct(:,1)),[length(bAct(1,:)) 1]);
               h=errorbar(handle,Xaxis',bAct,seAct,'.');%strcat('.',FitResult.colors{mod(i-1,length(FitResult.colors))+1})); 
               set(h, 'LineStyle', 'none', 'Marker', '.');%,...
               
               for n=1:length(h)
                   currColor = get(h(n),'Color');
                   set(h(n),'MarkerEdgeColor',currColor,'MarkerFaceColor',currColor); 

               end
            end

            hold on;
           
            
            if(plotSignificance==1)
               v=axis;
               vdiff = .8*v(4);

               for i=1:length(fitNum)
                  plot(handle,find(sigIndex(:,i)==1),vdiff*ones(length(find(sigIndex(:,i)==1)),1)-i*.1,strcat('*',FitResult.colors{mod(i-1,length(FitResult.colors))+1})); hold on;
               end
            end
           ylabel('GLM Fit Coefficients','Interpreter','none');
           xtickLabels = fitObj.getPlotParams.xLabels(subIndex);
           xticks = 1:(length(xtickLabels));
           
           set(handle,'xtick',xticks,'xtickLabel',xtickLabels,'FontSize',6);
%            axis tight;
           if(max(fitObj.numCoeffs)>=1)
            xticklabel_rotate([],90,[],'Fontsize',10);
           end
%            hT=rotateticklabel(gca,-90);
           h_legend=legend(handle,fitObj.lambda.dataLabels(fitNum),'Location','NorthEast');
           set(h_legend,'FontSize',14)
            pos = get(h_legend,'position');
            set(h_legend, 'position',[pos(1)+.05 pos(2) pos(3:4)]);
            
%            axis tight;
            title({'GLM Coefficients with 95% CIs (* p<0.05)'},'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Arial');
            set(gca,'FontName', 'Arial' );
            set(gca, ...
              'TickLength'  , [.02 .02] , ...
              'YGrid'       , 'on'      , ...
              'LineWidth'   , 1         );
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
           

        end
        function h=plotResults(fitObj)
            % plotResults(fitObj)
            % Generates KS plot, auto-correlation function of the inverse
            % gaussian transformed rescaled ISIs, the sequential
            % correlation coefficient between neigboring pairs of the
            % rescaled ISIs (zj vs. zj-1), the GLM regression coefficients,
            % and the Point Process Residual.
                scrsz = get(0,'ScreenSize');
                h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
                
                subplot(2,4,[1 2]); fitObj.KSPlot; %make the plot
                ht=text(.45, .95,strcat('Neuron:',num2str(fitObj.neuronNumber)));
                set(ht,'FontName', 'Arial','FontWeight','bold','FontSize',10);

                subplot(2,4,3); fitObj.plotInvGausTrans;
                subplot(2,4,4); fitObj.plotSeqCorr;
                subplot(2,4,[7 8]); fitObj.plotResidual;
                subplot(2,4,[5 6]); fitObj.plotCoeffs;
        end
        function handle = KSPlot(fitObj,fitNum)
            % handle = KSPlot(fitObj)
            % computes the K-S plot for each of the the candidate rate
            % functions in this FitResult object. These candidate rate
            % functions are numbered according to the order in which they
            % were added to the FitResult.
            if(nargin<2)
                fitNum=1:fitObj.numResults;
            end
            h=gcf; 
            %h=[];
            figure(h); 
        %     size(xAxis) 
        %     size(KSSorted)
            N = length(fitObj.KSStats.KSSorted);
            xaxis = fitObj.KSStats.xAxis(:,1);
            % Plot the CIs
            plot(xaxis,xaxis, 'k-.'); hold on;
            plot(xaxis, xaxis+1.36/sqrt(N), 'r','Linewidth',1); 
            plot(xaxis,xaxis-1.36/sqrt(N), 'r','Linewidth',1 );
            handle=plot(fitObj.KSStats.xAxis(:,fitNum),fitObj.KSStats.KSSorted(:,fitNum),'Linewidth',2);
            
            %set(gca,'xtick',[],'ytick',[],'ztick', [])
            axis( [0 1 0 1] );
%             dataLabels = cell(1,fitObj.lambda.dimension);
%             for i=1:fitObj.lambda.dimension
                dataLabels = fitObj.lambda.dataLabels(fitNum);
%             end
            h_legend=legend(handle,dataLabels,'Location','SouthEast');
            set(h_legend,'FontSize',14)
            hx=xlabel('Ideal Uniform CDF');
            hy=ylabel('Empirical CDF');
            title({'KS Plot of Rescaled ISIs'; 'with 95% Confidence Intervals'},'FontWeight','bold','FontSize',11,'FontName','Arial');
            set([hx, hy],'FontName', 'Arial','FontWeight','bold','FontSize',12);

            set(gca, ...
              'TickLength'  , [.02 .02] , ...
              'YTick'       , 0:.2:1, ...
              'XTick'       , 0:.2:1, ...
              'LineWidth'   , 1         );
            end

        function structure = toStructure(fitObj)
            % structure =  toStructure(fitObj)
            % Converts FitResult object to a matlab structure than can then
            % be saved. The structure is compatible with the FitResult
            % static method FitResult.fromStructure(structure) that returns
            % the object corresponding structure passed in.
            fnames = fieldnames(fitObj);
            
            for i=1:length(fnames)
                
                currObj = fitObj.(fnames{i});
                if(strcmp(fnames{i},'histObjects')||strcmp(fnames{i},'ensHistObjects'))
                    for j=1:fitObj.numResults
                        tempObj = fitObj.(fnames{i}){j};
                        if(~isempty(tempObj))
                            structure.(fnames{i}){j} = tempObj.toStructure;
                        else
                            structure.(fnames{i}){j} = tempObj;
                        end
                    end
                elseif(strcmp(fnames{i},'invGausStats'))
                    tempNames = fieldnames(fitObj.(fnames{i}));
                    for j=1:length(tempNames)
                       tempObj = currObj.(tempNames{j});
                       if(~isempty(tempObj))
                        structure.(fnames{i}).(tempNames{j})=  tempObj.dataToStructure;
                       else
                         structure.(fnames{i}).(tempNames{j})=  tempObj;  
                       end
                        
                    end
                    
                else
                
                    if(isa(currObj,'double')||isa(currObj,'cell'))
                        structure.(fnames{i}) = currObj;
                    elseif(isa(currObj,'Covariate') ||isa(currObj,'ConfigColl')||isa(currObj,'nspikeTrain'))
                        structure.(fnames{i}) = currObj.toStructure;
                    elseif(isa(currObj,'SignalObj'))
                        structure.(fnames{i})  = currObj.dataToStructure;
                    elseif(isa(currObj,'struct'))
                        structure.(fnames{i}) = currObj;
                    end
                end
            end
            
        end
        

        function handle = plotSeqCorr(fitObj)
            % handle = plotSeqCorr(fitObj)
            % plot zj+1 against zj
            
            %colors = {'.b','.g','.r','.c','.m','.y','.k'};
            rho=zeros(1,fitObj.numResults);
            pval=zeros(1,fitObj.numResults);
            dataLabels = fitObj.lambda.dataLabels;
            for i=1:fitObj.numResults
               handle = plot(fitObj.U(1:end-1,i),fitObj.U(2:end,i),strcat('.',Analysis.colors{mod(i-1,length(Analysis.colors))+1})); hold on;
               [rhoTemp,p]= corrcoef(fitObj.U(1:end-1,i),fitObj.U(2:end,i));%handle=scatterhist(fitResults.Z(1:end-1,i),fitResults.Z(2:end,i))
               
               [~,columns]=size(rhoTemp);
                if(columns>1)
                    rho(i) = rhoTemp(1,2);
                    pval(i)= p(1,2);
                else
                    rho(i) = rhoTemp;
                    pval(i)= p;
                end
               dataLabels{i} = strcat(dataLabels{i},', \rho=',num2str(rho(i),'%0.2g'),' (p=',num2str(pval(i),'%0.2g'),')');
                %get(h,'AlphaData');
                %set(h,'FaceAlpha',0.2,'EdgeAlpha',0.8,'EdgeColor',color{i});
            end
            
           
            h_legend=legend(dataLabels,'Location','NorthEast');     
            set(h_legend,'FontSize',14)
            pos = get(h_legend,'position');
            set(h_legend, 'position',[pos(1)+.05 pos(2) pos(3:4)]);
            
            hy=ylabel('u_{j+1}'); hx=xlabel('u_j');
            set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');

            axis([0 1 0 1]);
            title({'Sequential Correlation of'; 'Rescaled ISIs'},'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Arial');
        
    
            set(gca, ...
              'TickLength'  , [.02 .02] , ...
              'YTick'       , 0:.25:1, ...
              'XTick'       , 0:.25:1, ...
              'LineWidth'   , 1         );
        

            
        end
        function handle = plotInvGausTrans(fitObj)
            % handle = plotInvGausTrans(fitObj)
            % Plots the Auto-correlation function of the X_j's where:
            % Z_j: rescaled ISI from the Time Rescaling Theorem. 
            %      Exponential Rate 1 under true conditional intensity
            %      function.
            % U_j: 1-exp(-Z_j). Uniform on the interval [0,1) if Z_j's
            %      are exponential rate 1.
            %
            % X_j: norminv(U_j,0,1). Gaussian mean 0, stdev 1 if U_j's are
            %      U([0,1))
            %
            
            %[rows,colm] = size(fitObj.X);
            %index=find(fitObj.invGausStats.lags==1);
            %lags=fitObj.invGausStats.lags;
            [fitObj.X,rhoSig,confBoundSig] = Analysis.computeInvGausTrans(fitObj.Z);
%             rhoSig=fitObj.invGausStats.rhoSig;
            n=length(fitObj.X);
%             confBoundSig = fitObj.invGausStats.confBoundSig;
            handle=[];
%              for i=1:colm
%                     %i
%                     htemp=plot(lags',rho(:,i),strcat('.',FitResults.colors{mod(i-1,length(Analysis.colors))+1}));
%                     handle=[handle,htemp];
%                     hold on; 
%                     %labelArray{i} = ['Fit ' num2str(i)];
%              end
            
            rhoSig.plot;
            h_legend=legend(fitObj.lambda.dataLabels,'Location','NorthEast');            
            set(h_legend,'FontSize',14)
            pos = get(h_legend,'position');
            set(h_legend, 'position',[pos(1)+.05 pos(2) pos(3:4)]);
            %legend(h,labelArray); 
            hold on; confBoundSig.plot;
            title({'Autocorrelation Function';'of Rescaled ISIs'; 'with 95% CIs'},'FontWeight','bold',...
            'FontSize',11,...
            'FontName','Arial');
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
            set(gca, ...
              'TickLength'  , [.02 .02] , ...
              'LineWidth'   , 1         );
            v=axis;
            maxY = max(abs(v(3:4)))*(1.1); %add 10%
            axis([v(1:2) -maxY maxY]);
        end
        function handle = plotResidual(fitObj)
            % handle = plotResidual(fitObj)
            % Plots the Point Process Residual
            handle=fitObj.Residual.plot;
            legend off;
            h_legend=legend(fitObj.lambda.dataLabels,'Location','NorthEast'); 
            set(h_legend,'FontSize',14)
            pos = get(h_legend,'position');
%             set(h_legend, 'position',[.91 .41 pos(3:4)]);
            set(h_legend, 'position',[pos(1)+.05 pos(2) pos(3:4)]);
            title('Point Process Residual','FontWeight','bold',...
            'FontSize',11,...
            'FontName','Arial');
            xlabel('time [s]','Interpreter','none');
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
            v=axis;
            maxY = max(abs(v(3:4)))*(1.1); %add 10%
            axis([v(1:2) -maxY maxY]);
        end
        
        
        function setKSStats(fitObj, Z, U, xAxis, KSSorted, ks_stat)
            % setKSStats(fitObj, Z, xAxis, KSSorted, ks_stat)
            % Allows KS statistics to be set after object creation
            % Z: Rescaled ISIs from the Time Rescaling Theorem
            % xAxis: xAxis of the KS plot
            % KSSorted: the sorted values of Uj=1-exp(-zj)
            % ks_stat: the maximum deviation from the 45 degree line for
            % all of the fits.
            % 
            fitObj.Z        =Z;
            fitObj.U        =U;
            fitObj.KSStats.xAxis    =xAxis;
            fitObj.KSStats.KSSorted =KSSorted;
            
            for i=1:size(xAxis,2);
                [differentDists(i),pVal(i),ks_stat(i)]=kstest2(fitObj.KSStats.xAxis(:,i) ,fitObj.KSStats.KSSorted(:,i));
            end
            
            fitObj.KSStats.ks_stat  =ks_stat;
%             N = length(fitObj.KSStats.KSSorted);
%             fitObj.KSStats.withinConfInt = ks_stat<1.36/sqrt(N);
            fitObj.KSStats.withinConfInt = ~differentDists;
            fitObj.KSStats.pValue = pVal;
        end
        function setInvGausStats(fitObj, X,rhoSig,confBoundSig)
            % setInvGausStats(fitObj,X,rhoSig,confBoundSig)
            % Sets the inverse gaussian transformed rescaled ISIs and the
            % confidence bounds after the object has been created.
            %fitObj.U=U;
            fitObj.X=X;
            fitObj.invGausStats.rhoSig=rhoSig;
            fitObj.invGausStats.confBoundSig=confBoundSig;
        end        
        function setFitResidual(fitObj,M)
           % setFitResidual(fitObj,M). 
           % Adds the point process residual to the FitResult object
           fitObj.Residual = M; 
        end
        
        function [paramVals, paramSE, paramSigIndex] = getParam(fitObj,paramNames,fitNum)
           % output is a matrix of length equal to the total number of
           % paramNames
           % and one column for each fit
             
            if(nargin<3)
                fitNum = 1:fitObj.numResults;
            end
           if(isempty(fitObj.plotParams))
               fitObj.computePlotParams;
           end
           paramVals = zeros(length(paramNames),length(fitNum));
           
           if(nargout>1)
            paramSE = zeros(length(paramNames),length(fitNum));
           end
           
           if(nargout>2)
            paramSigIndex = zeros(length(paramNames),length(fitNum));
           end
           
           for i=1:length(paramNames)
                paramIndex=find(strcmp(paramNames(i),fitObj.uniqueCovLabels));
                paramVals(i,:) = fitObj.plotParams.bAct(paramIndex,fitNum);
                if(nargout>1)
                    paramSE(i,:) = fitObj.plotParams.seAct(paramIndex,fitNum);
                end
                if(nargout>2)
                    paramSigIndex(i,:) = fitObj.plotParams.sigIndex(paramIndex,fitNum);
                end
           end
            
        end
        
            
    end
    
    methods (Static)
        
        function fitObj = fromStructure(structure)
            % fitObj = fromStructure(structure)
            % Returns a FitResult object from the input structure.
            % This is used to be able to save and restore FitResult
            if(isa(structure,'struct'))
                if(isa(structure.neuralSpikeTrain,'cell'))
                    spikeObj = cell(1,length(structure.neuralSpikeTrain));
                   for k=1:length(structure.neuralSpikeTrain)
                       spikeObj{k} = nspikeTrain.fromStructure(structure.neuralSpikeTrain{k}); 
                       
                   end
                else
                    spikeObj=nspikeTrain.fromStructure(structure.neuralSpikeTrain); 
                end
                lambda=Covariate.fromStructure(structure.lambda);
                rhoSig=SignalObj.signalFromStruct(structure.invGausStats.rhoSig);
                confBoundSig = SignalObj.signalFromStruct(structure.invGausStats.confBoundSig);
                M = Covariate.fromStructure(structure.Residual);
                for i=1:structure.numResults
                    histObjects{i} = History.fromStructure(structure.histObjects{i});
                    ensHistObject{i} = History.fromStructure(structure.ensHistObjects{i});
                end
                configColl = ConfigColl.fromStructure(structure.configs);
                fitObj=FitResult(spikeObj,structure.covLabels,structure.numHist,histObjects,ensHistObject,lambda,structure.b, structure.dev, structure.stats,structure.AIC,structure.BIC,configColl,structure.XvalData,structure.XvalTime,structure.fitType);
                fitObj.setKSStats(structure.Z,structure.U, structure.KSStats.xAxis, structure.KSStats.KSSorted, structure.KSStats.ks_stat);
                fitObj.setInvGausStats(structure.X,rhoSig,confBoundSig);
                fitObj.setFitResidual(M);
                fitObj.setNeuronName(structure.neuronNumber);
                
            elseif(isa(structure,'cell')) %cell array of FitResult objects
                fitObj = cell(size(structure));
               for i=1:length(structure)
                  fitObj{i} = FitResult.fromStructure(structure{i}); 
               end
            end
            
        end
        
        function structCell = CellArrayToStructure(fitResObjCell)
            % structCell = CellArrayToStructure(fitResObjCell)
            % For every FitResult structure in the a cell array, it calls 
            % FitResults.fromStructure
            
           if(isa(fitResObjCell,'FitResult'))
               structCell = fitResObjCell.toStructure;
           elseif(isa(fitResObjCell,'cell')&&~isempty(fitResObjCell))
               if(isa(fitResObjCell{1},'FitResult'))
                   structCell = cell(size(fitResObjCell));
                  for i=1:length(fitResObjCell)
                     structCell{i} = fitResObjCell{i}.toStructure;
                  end
               end
           else
            structCell={};
           end
            
        end
        
    end
    
end

%Helper functions
function hText = xticklabel_rotate(XTick,rot,varargin)
    %hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
    %
    % Syntax: xticklabel_rotate
    %
    % Input:    
    % {opt}     XTick       - vector array of XTick positions & values (numeric) 
    %                           uses current XTick values or XTickLabel cell array by
    %                           default (if empty) 
    % {opt}     rot         - angle of rotation in degrees, 90° by default
    % {opt}     XTickLabel  - cell array of label strings
    % {opt}     [var]       - "Property-value" pairs passed to text generator
    %                           ex: 'interpreter','none'
    %                               'Color','m','Fontweight','bold'
    %
    % Output:   hText       - handle vector to text labels
    %
    % Example 1:  Rotate existing XTickLabels at their current position by 90°
    %    xticklabel_rotate
    %
    % Example 2:  Rotate existing XTickLabels at their current position by 45° and change
    % font size
    %    xticklabel_rotate([],45,[],'Fontsize',14)
    %
    % Example 3:  Set the positions of the XTicks and rotate them 90°
    %    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
    %    xticklabel_rotate([1960:2:2004]);
    %
    % Example 4:  Use text labels at XTick positions rotated 45° without tex interpreter
    %    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
    %
    % Example 5:  Use text labels rotated 90° at current positions
    %    xticklabel_rotate([],90,NameFields);
    %
    % Note : you can not re-run xticklabel_rotate on the same graph. 
    %
    % 


    % This is a modified version of xticklabel_rotate90 by Denis Gilbert
    % Modifications include Text labels (in the form of cell array)
    %                       Arbitrary angle rotation
    %                       Output of text handles
    %                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size 
    %                           and keep text on plot
    %                           (handles small window resizing after, but not well due to proportional placement with 
    %                           fixed font size. To fix this would require a serious resize function)
    %                       Uses current XTick by default
    %                       Uses current XTickLabel is different from XTick values (meaning has been already defined)

    % Brian FG Katz
    % bfgkatz@hotmail.com
    % 23-05-03
    % Modified 03-11-06 after user comment
    %	Allow for exisiting XTickLabel cell array

    % Other m-files required: cell2mat
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: xticklabel_rotate90, TEXT,  SET

    % Based on xticklabel_rotate90
    %   Author: Denis Gilbert, Ph.D., physical oceanography
    %   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
    %   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
    %   February 1998; Last revision: 24-Mar-2003

    % check to see if xticklabel_rotate has already been here (no other reason for this to happen)
    if isempty(get(gca,'XTickLabel')),
        error('xticklabel_rotate : can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased')  ;
    end

    % if no XTickLabel AND no XTick are defined use the current XTickLabel
    %if nargin < 3 & (~exist('XTick') | isempty(XTick)),
    if (nargin < 3 || isempty(varargin{1})) && (~exist('XTick') || isempty(XTick)),
        xTickLabels = get(gca,'XTickLabel')  ; % use current XTickLabel
        if ~iscell(xTickLabels)
            % remove trailing spaces if exist (typical with auto generated XTickLabel)
            temp1 = num2cell(xTickLabels,2)         ;
            for loop = 1:length(temp1),
                temp1{loop} = deblank(temp1{loop})  ;
            end
            xTickLabels = temp1                     ;
        end
    varargin = varargin(2:length(varargin));	
    end

    % if no XTick is defined use the current XTick
    if (~exist('XTick') | isempty(XTick)),
        XTick = get(gca,'XTick')        ; % use current XTick 
    end

    %Make XTick a column vector
    XTick = XTick(:);

    if ~exist('xTickLabels'),
        % Define the xtickLabels 
        % If XtickLabel is passed as a cell array then use the text
        if (length(varargin)>0) & (iscell(varargin{1})),
            xTickLabels = varargin{1};
            varargin = varargin(2:length(varargin));
        else
            xTickLabels = num2str(XTick);
        end
    end    

    if length(XTick) ~= length(xTickLabels),
        error('xticklabel_rotate : must have same number of elements in "XTick" and "XTickLabel"')  ;
    end

    %Set the Xtick locations and set XTicklabel to an empty string
    set(gca,'XTick',XTick,'XTickLabel','')

    if nargin < 2,
        rot = 90 ;
    end

    % Determine the location of the labels based on the position
    % of the xlabel
    hxLabel = get(gca,'XLabel');  % Handle to xlabel
    xLabelString = get(hxLabel,'String');

    % if ~isempty(xLabelString)
    %    warning('You may need to manually reset the XLABEL vertical position')
    % end

    set(hxLabel,'Units','data');
    xLabelPosition = get(hxLabel,'Position');
    y = xLabelPosition(2);

    %CODE below was modified following suggestions from Urs Schwarz
    y=repmat(y,size(XTick,1),1);
    % retrieve current axis' fontsize
    fs = get(gca,'fontsize');

    % Place the new xTickLabels by creating TEXT objects
    hText = text(XTick, y, xTickLabels,'fontsize',fs);

    % Rotate the text objects by ROT degrees
    set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})

    % Adjust the size of the axis to accomodate for longest label (like if they are text ones)
    % This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
    % This approach keeps the right side of the graph at the same place 

    set(get(gca,'xlabel'),'units','data')           ;
        labxorigpos_data = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'ylabel'),'units','data')           ;
        labyorigpos_data = get(get(gca,'ylabel'),'position')  ;
    set(get(gca,'title'),'units','data')           ;
        labtorigpos_data = get(get(gca,'title'),'position')  ;

    set(gca,'units','pixel')                        ;
    set(hText,'units','pixel')                      ;
    set(get(gca,'xlabel'),'units','pixel')          ;
    set(get(gca,'ylabel'),'units','pixel')          ;

    origpos = get(gca,'position')                   ;
    
    % Modified by Iahn Cajigas
    % 3/4/2011 to allow for a single
    % xTickLabel to work properly
    temphText = get(hText,'extent');
    if(isa(temphText,'cell'))
        textsizes = cell2mat(temphText)       ;
    else
        textsizes = temphText;
    end
    
    longest =  max(textsizes(:,4))                  ;

    laborigext = get(get(gca,'xlabel'),'extent')    ;
    laborigpos = get(get(gca,'xlabel'),'position')  ;


    labyorigext = get(get(gca,'ylabel'),'extent')   ;
    labyorigpos = get(get(gca,'ylabel'),'position') ;
    leftlabdist = labyorigpos(1) + labyorigext(1)   ;

    % assume first entry is the farthest left
    leftpos = get(hText(1),'position')              ;
    leftext = get(hText(1),'extent')                ;
    leftdist = leftpos(1) + leftext(1)              ;
    if leftdist > 0,    leftdist = 0 ; end          % only correct for off screen problems

    botdist = origpos(2) + laborigpos(2)            ;
    newpos = [origpos(1)-leftdist longest+botdist origpos(3)+leftdist origpos(4)-longest+origpos(2)-botdist]  ;
    set(gca,'position',newpos)                      ;

    % readjust position of nex labels after resize of plot
    set(hText,'units','data')                       ;
    for loop= 1:length(hText),
        set(hText(loop),'position',[XTick(loop), y(loop)])  ;
    end


    % adjust position of xlabel and ylabel
    laborigpos = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0])   ;

    % switch to data coord and fix it all
    set(get(gca,'ylabel'),'units','data')                   ;
    set(get(gca,'ylabel'),'position',labyorigpos_data)      ;
    set(get(gca,'title'),'position',labtorigpos_data)       ;

    set(get(gca,'xlabel'),'units','data')                   ;
        labxorigpos_data_new = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)])   ;


    % Reset all units to normalized to allow future resizing
    set(get(gca,'xlabel'),'units','normalized')          ;
    set(get(gca,'ylabel'),'units','normalized')          ;
    set(get(gca,'title'),'units','normalized')          ;
    set(hText,'units','normalized')                      ;
    set(gca,'units','normalized')                        ;

    if nargout < 1,
        clear hText
    end

end
function [uniqueLabels, indexIntoOriginal, restoreIndex] = getUniqueLabels(covLabels)
% Given a set of covLabels, returns a subset of labels that are unique
            offset = 0;
            for i=1:length(covLabels)
                currLabels = covLabels{i};                
                allLabels((1:length(currLabels))+offset) = currLabels;
                offset=length(allLabels);
            end
            [uniqueLabels, indexIntoOriginal, restoreIndex] = unique(allLabels);
end
        

   