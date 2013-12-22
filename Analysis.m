classdef Analysis
% ANALYSIS Collection of functions (static methods) used for GLM analysis
% of point process data.
% <a href="matlab: methods('Analysis')">methods</a>
% <a href="matlab:web('AnalysisExamples.html', '-helpbrowser')">Analysis Examples</a> 
%
% see also <a href="matlab:help('Trial')">Trial</a>, <a
% href="matlab:help('CovColl')">CovColl</a>, <a
% href="matlab:help('nstColl')">nstColl</a>,<a
% href="matlab:help('History')">History</a>
%
% Reference page in Help browser
% <a href="matlab: doc('Analysis')">doc Analysis</a>

%%
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



properties (Constant)
    colors = {'b','g','r','c','m','y','k'};
end

    methods (Static)
        function fitResults =RunAnalysisForNeuron(tObj,neuronNumber,configColl,makePlot,Algorithm,DTCorrection,batchMode)
            % fitResults =RunAnalysisForNeuron(tObj,neuronNumber,configColl,makePlot,Algorithm)
            % tObj: Trial to be analyzed
            % neuronNumber: number of the neuron to be analyzed. Can be a  
            %               vector to specify multiple neurons to be analyzed.
            %               If more than one neuron specified, then
            %               fitResults is a cell array of fitResult
            %               objects. fitResults{i} will contain the
            %               fitResults object for neuronNum(i).
            % configColl: ConfigColl object containing the different
            %             configurations (description of the the types of fits, eg. covariates) that correspond to each fit.
            % makePlot: Set to 1 to show a summary plot for this neuron. If performing multiple neuron analysis (eg. via RunAnalysisForAllNeurons) set ths parameter to zero to avoid screen clutter.
            % Algorithm: Either 'GLM' or 'BNLRCG'. Default is 'GLM'
            %            GLM - Standard GLM regression from matlab.
            %            BNLRCG - faster Truncated, L-2 Regularized,
            %            Binomial Logistic Regression with Conjugate
            %            Gradient Solver by Demba Ba (demba@mit.edu).
            % DTCorrection: 0 for no DT Correction of KS plot, 1 is the
            %               default.
            %
            % batchMode: when batchMode=1 neurons with same name are fit at once rather than separetely
        
            if(nargin<7 || isempty(batchMode))
                batchMode = 0; %default treat each spike train separately
            end
            
            if(nargin<6 || isempty(DTCorrection))
                DTCorrection =1;
            end
            if(nargin<5 || isempty(Algorithm))
                Algorithm = 'GLM';
            end
            if(nargin<4 || isempty(makePlot))
                makePlot=1;
            end
            numNeurons = length(neuronNumber);
            labels=cell(numNeurons,1);
            lambda=cell(numNeurons,1);
            b     =cell(numNeurons,1);
            dev   =zeros(numNeurons,1);
            numHist=cell(numNeurons,1);
            stats =cell(numNeurons,1);
            histObj =cell(numNeurons,1);
            ensHistObj=cell(numNeurons,1);
            AIC   =zeros(numNeurons,1);
            BIC   =zeros(numNeurons,1);
            windowSize = .01; % 1/tObj.sampleRate;% for Residual Computation;
            spikeTraining = cell(1,numNeurons);
            XvalData =cell(numNeurons,1);
            XvalTime =cell(numNeurons,1);
            spikeValidation = cell(1,numNeurons);
            %% Fit Using Training Data
            if(diff(tObj.validationWindow)~=0)
                tObj.setTrialTimesFor('training');
            end
            if(batchMode==1)
               display('Running in batch mode: neurons with same name are fit simultaneously'); 
            end
            for i=1:configColl.numConfigs
                configColl.setConfig(tObj,i);
%                 fprintf(strcat('Analyzing Configuration #',num2str(i)));
                
                if(batchMode==0)
                    fprintf(strcat('Analyzing Configuration #',num2str(i),': Neuron #'));
                    for j=1:numNeurons
%                         fprintf(strcat('Analyzing Configuration #',num2str(i),': Neuron #',num2str(neuronNumber(j))));
                        if(j==1)
                            fprintf('%d',neuronNumber(j));
                        else
                            fprintf(',%d',neuronNumber(j));
                        end
                        %clear tempLabels;
                        %tObj.setCurrentNeuron(neuronNumber);
                        otherLabels  = tObj.getLabelsFromMask(neuronNumber(j));
    %                     labels{j}{i}  = horzcat('Baseline',otherLabels); % Labels change depending on presence/absense of History or ensCovHist
                        labels{j}{i}  = otherLabels; % Labels change depending on presence/absense of History or ensCovHist
                        numHist{j}{i} = tObj.getNumHist;
                        histObj{j}{i} = tObj.history;
                        ensHistObj{j}{i} = tObj.ensCovHist;
                        [lambdaTemp, bTemp, devTemp, statsTemp,AICTemp,BICTemp,distribTemp] = Analysis.GLMFit(tObj,neuronNumber(j),i,Algorithm);
                        lambda{j}{i} = lambdaTemp; b{j}{i} = bTemp; stats{j}{i} = statsTemp;
                        dev(j,i) = devTemp;  AIC(j,i)= AICTemp; BIC(j,i)= BICTemp;
                        distrib{j}{i} =distribTemp;
                        spikeTraining{j} = tObj.nspikeColl.getNST(neuronNumber(j));%.nstCopy;
                        spikeTraining{j}.setName(num2str(neuronNumber(j)));

                        %% Collect the validation Data
                        if(diff(tObj.validationWindow)~=0)
                              tObj.setTrialTimesFor('validation');
                              XvalData{j}{i}=tObj.getDesignMatrix(neuronNumber(j));
                              XvalTime{j}{i}=tObj.covarColl.getCov(1).time;
                              spikeValidation{j} = tObj.nspikeColl.getNST(neuronNumber(j));%.nstCopy;
                              spikeTraining{j}.setName(num2str(neuronNumber(j)));
                              tObj.setTrialTimesFor('training')
                        end
                       
                    end
                elseif(batchMode==1)
                    neuronNames=neuronNumber; % This is an index of names in the batchMode case
                    fprintf(strcat('Analyzing Configuration #',num2str(i),': Neuron #'));
                    for j=1:numNeurons
                        
%                         if(isa(neuronNames,'cell'))
%                             fprintf(strcat('Analyzing Configuration #',num2str(i),': Neuron #'));
%                             display(strcat('Analyzing Configuration #',num2str(i),': Neuron #',neuronNames{j}));
%                         elseif(isa(neuronNames,'char'))
%                             display(strcat('Analyzing Configuration #',num2str(i),': Neuron #',neuronNames));
%                         elseif(isa(neuronNames,'double'))
%                             display(strcat('Analyzing Configuration #',num2str(i),': Neuron #',num2str(neuronNames)));
%                         end
                        if(isa(neuronNames,'cell'))
                            if(j==1)
                                fprintf('%s',neuronNames{j});
                            else
                                fprintf(',%s',neuronNames{j});
                            end
                        elseif(isa(neuronNames,'char'))
                            if(j==1)
                                fprintf('%s',neuronNames);
                            else
                                fprintf(',%s',neuronNames);
                            end
                        elseif(isa(neuronNames,'double'))
                            if(j==1)
                                fprintf('%d',neuronNames);
                            else
                                fprintf('%,d',neuronNames);
                            end
                        end

                        %clear tempLabels;
                        %tObj.setCurrentNeuron(neuronNumber);
                        otherLabels  = tObj.getLabelsFromMask(neuronNumber(j));
    %                     labels{j}{i}  = horzcat('Baseline',otherLabels); % Labels change depending on presence/absense of History or ensCovHist
                        labels{j}{i}  = otherLabels; % Labels change depending on presence/absense of History or ensCovHist
                        numHist{j}{i} = tObj.getNumHist;
                        histObj{j}{i} = tObj.history;
                        ensHistObj{j}{i} = tObj.ensCovHist;
                         if(isa(neuronNames,'cell'))
                             [lambdaTemp, bTemp, devTemp, statsTemp,AICTemp,BICTemp,distribTemp] = Analysis.GLMFit(tObj,neuronNumber{j},i,Algorithm);
                         elseif(isa(neuronNames,'char'))
                             [lambdaTemp, bTemp, devTemp, statsTemp,AICTemp,BICTemp,distribTemp] = Analysis.GLMFit(tObj,neuronNumber,i,Algorithm);
                         else
                            [lambdaTemp, bTemp, devTemp, statsTemp,AICTemp,BICTemp,distribTemp] = Analysis.GLMFit(tObj,neuronNumber(j),i,Algorithm);
                         end
                       
                        lambda{j}{i} = lambdaTemp; b{j}{i} = bTemp; stats{j}{i} = statsTemp;
                        dev(j,i) = devTemp;  AIC(j,i)= AICTemp; BIC(j,i)= BICTemp;
                        distrib{j}{i} =distribTemp;
                         if(isa(neuronNames,'cell'))
                            currSpikes=tObj.nspikeColl.getNST(tObj.getNeuronIndFromName(neuronNames{j}));
                         elseif(isa(neuronNames,'char'))
                            currSpikes=tObj.nspikeColl.getNST(tObj.getNeuronIndFromName(neuronNames));
                         else
                            currSpikes=tObj.nspikeColl.getNST(neuronNames(j)); 
                         end
                        
                        
                        for n=1:length(currSpikes)
                            if(isa(currSpikes,'cell'))
                                 currSpikes{n} = currSpikes{n}.nstCopy;
                                 if(isa(neuronNames,'cell'))
                                    currSpikes{n}.setName(neuronNames{j});
                                 elseif(isa(neuronNames,'char'))
                                    currSpikes{n}.setName(neuronNames);
                                 else
                                    currSpikes{n}.setName(neuronNames(j));
                                 end
                                        
                            else
                                currSpikes = currSpikes.nstCopy;
                              
                                if(isa(neuronNames,'cell'))
                                    currSpikes.setName(neuronNames{j});
                                elseif(isa(neuronNames,'char'))
                                    currSpikes.setName(neuronNames);
                                else
                                    currSpikes.setName(neuronNames(j));
                                end
                                 
                            end
                        end
                            
                        spikeTraining{j} = currSpikes;

                        %% Collect the validation Data
                        if(diff(tObj.validationWindow)~=0)
                            tObj.setTrialTimesFor('validation');
                            tempIndices=tObj.getNeuronIndFromName(neuronNames{j});
                            currSpikes=tObj.nspikeColl.getNST(tempIndices);
                            tempX = [];
                            tempTime=[];
                            for n=1:length(tempIndices)
                                currSpikes{n} = currSpikes{n}.nstCopy;
                                currSpikes{n}.setName(neuronNames{j});
                                if(n==1)
                                    tempX =tObj.getDesignMatrix(tempIndices(n)); 
                                    tempTime =tObj.covarColl.getCov(1).time;
                                else
                                    tempX = [tempX; tObj.getDesignMatrix(tempIndices(n))];
                                    offset = max(tempTime)+1/tObj.sampeRate;
                                    tempTime = [tempTime;(tObj.covarColl.getCov(1).time+offset)];
                                end
                            end
                            spikeValidation{j} = currSpikes;
                            XvalData{j}{i}=tempX;
                            XvalTime{j}{i}=tempTime;
                            
                            tObj.setTrialTimesFor('training')
                        end
                    end
                    
                    
                end
                fprintf('\n');
            end
            
            
%             %% Collect the validation Data
% 
%             if(diff(tObj.validationWindow)~=0)
%                 tObj.setTrialTimesFor('validation');
%                 for i=1:configColl.numConfigs
%                     configColl.setConfig(tObj,i);
%                     for j=1:numNeurons
%                         XvalData{j,i}=tObj.getDesignMatrix(neuronNumber(j));
%                         XvalTime{j,i}=tObj.covarColl.getCov(1).time;
%                         spikeValidation{j} = tObj.nspikeColl.getNST(neuronNumber(j)).nstCopy;
%                         spikeTraining{j}.setName(num2str(neuronNumber(j)));
%                     end
%                 end
%                 
%                 %tObj.setTrialTimesFor('training');
%             end
%             
            
            %% Store the results

            fitResults =cell(length(neuronNumber),1);
            for j=1:numNeurons
                fitResults{j}=FitResult(spikeTraining{j},labels{j},numHist{j},histObj{j},ensHistObj{j},lambda{j},b{j}, dev(j,:), stats{j},AIC(j,:),BIC(j,:),configColl,XvalData{j},XvalTime{j},distrib{j});
                if(diff(tObj.validationWindow)~=0)
                    tObj.setTrialTimesFor('validation');
                    lambdaValidation = fitResults{j}.computeValLambda;
                    ValResults = FitResult(spikeValidation{j},labels{j},numHist{j},histObj{j},ensHistObj{j},lambdaValidation,b{j}, dev(j,:), stats{j},AIC(j,:),BIC(j,:),configColl,XvalData{j},XvalTime{j},distrib);
                    fitResults{j}.validation = ValResults; %validation field is actually another fitResults object but with the validation data
                end
                %% Process the results and compute further parameters
                if(makePlot==1)
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[scrsz(3)*.1 scrsz(4)*.1 scrsz(3)*.8 scrsz(4)*.8]);
                    subplot(2,4,[1 2]);     Analysis.KSPlot(fitResults{j},DTCorrection,makePlot); %make the plot 
                    hold on; text(.45, .95,strcat('Neuron:',num2str(neuronNumber(j))));
                    subplot(2,4,3);         Analysis.plotInvGausTrans(fitResults{j},makePlot);
                    subplot(2,4,4);         Analysis.plotSeqCorr(fitResults{j});
                    subplot(2,4,[7 8]);     Analysis.plotFitResidual(fitResults{j},windowSize,makePlot);
                    subplot(2,4,[5 6]);     Analysis.plotCoeffs(fitResults{j});
                else
                    Analysis.KSPlot(fitResults{j},DTCorrection,makePlot);
                    Analysis.plotInvGausTrans(fitResults{j},makePlot);
                    Analysis.plotFitResidual(fitResults{j},windowSize,makePlot);
                    %fitResults.computePlotParams;
                end
            end
            if(length(neuronNumber)==1)
                fitResults = fitResults{1};
            end

        end
        function fitResults = RunAnalysisForAllNeurons(tObj,configs,makePlot,Algorithm,DTCorrection,batchMode)
            % fitResults = RunAnalysisForAllNeurons(tObj,configs,makePlot,Algorithm)
            % Runs the fits specifed by configs (a ConfigColl object) on
            % all the neurons that are unmasked in the trial tObj. 
            % tObj - trial to be analyzed
            % configs - ConfigColl object specifying the types of fits to
            %           be performed.
            % makePlot - Set to 1 to generate a summary plot for each
            %            neuron.
            % Algorithm: Either 'GLM' or 'BNLRCG'. Default is 'GLM'
            %            GLM - Standard GLM regression from matlab.
            %            BNLRCG - faster Truncated, L-2 Regularized,
            %            Binomial Logistic Regression with Conjugate
            %            Gradient Solver by Demba Ba (demba@mit.edu).
            % DTCorrection: 0 for no DT Correction of KS plot, 1 is the
            %               default.
            % batchMode: when batchMode=1 neurons with same name are fit at once rather than separetely
        
            if(nargin<6 || isempty(batchMode))
                batchMode = 0; %default treat each spike train separately
            end
            if(nargin<5 || isempty(DTCorrection))
               DTCorrection =1; 
            end
            
            if(nargin<4 || isempty(Algorithm))
                Algorithm = 'GLM';
            end
            if(nargin<3 || isempty(makePlot))
                makePlot=1; %default to plotting results
            end
           
            
            if(batchMode==0)
                neuronIndex=tObj.getNeuronIndFromMask;
            elseif(batchMode==1)
                neuronIndex=tObj.getUniqueNeuronNames;
            end
%             numLoops = floor(length(neuronIndex)/4);
%             loopArray = cell(1,numLoops);
%             for k=1:numLoops
%                 if(k==numLoops)
%                     loopArray{k} = neuronIndex((4*(k-1)+1):end);
%                 else
%                     loopArray{k} = neuronIndex((4*(k-1)+1):4*k);
%                 end
%             end
            
           % parfor i=1:length(neuronIndex)
                fitResults = Analysis.RunAnalysisForNeuron(tObj,neuronIndex,configs,makePlot,Algorithm,DTCorrection,batchMode);
            %end
            
        end
            
        
        function [lambda,b, dev, stats,AIC, BIC,distribution] = GLMFit(tObj,neuronNumber,lambdaIndex,Algorithm)
            % [lambda,b, dev, stats,AIC, BIC] = GLMFit(tObj,neuronNumber,lambdaIndex,Algorithm)
            % Given a trial, tObj, and a neuronNumber specifying a neuron
            % from the trial, extracts the design matrix X from the current
            % covariate masks, history, and ensemble history in the trial,
            % and the observation vector,Y, and performs the GLM regression
            % using the specified algorithm. lambdaIndex: is used to
            % labeling the returned lambda with the number of the
            % configuration that it corresponds to. 
            % Algorithm: Either 'GLM' or 'BNLRCG'. Default is 'GLM'
            %            GLM - Standard GLM regression from matlab.
            %            BNLRCG - faster Truncated, L-2 Regularized,
            %            Binomial Logistic Regression with Conjugate
            %            Gradient Solver by Demba Ba (demba@mit.edu).
            % Returns:
            % lambda - Covariate containing the resulting conditional
            %          intensity function evaluated with the design matrix data.
            % b      - the GLM regression coefficients. Constant term is
            %          first followed by the components in X.
            %
            % dev    - deviance for the this regression.
            % stats  - stats structure from the GLM regression
            %          (p-values,std dev, etc.)
            % AIC    - Akaike's information criteria for this regression.
            % BIC    - Bayes Information Criteria for this regression.

            if(nargin<4)
              Algorithm='GLM';
            end
            if(isa(neuronNumber,'double'))
                binaryRep=tObj.nspikeColl.getNST(neuronNumber).isSigRepBinary;
                indices=neuronNumber;
            elseif(isa(neuronNumber,'char'))
                indices=tObj.getNeuronIndFromName(neuronNumber);
                binRep=zeros(size(indices));
                for i=1:length(indices)
                   binRep(i)=tObj.nspikeColl.getNST(indices(i)).isSigRepBinary;
                end
                binaryRep=prod(binRep);
            elseif(isa(neuronNumber,'cell'))
                indices=tObj.getNeuronIndFromName(neuronNumber{1});
                binRep=zeros(size(indices));
                for i=1:length(indices)
                   binRep(i)=tObj.nspikeColl.getNST(indices(i)).isSigRepBinary;
                end
                binaryRep=prod(binRep);
            end
            
            if(strcmp(Algorithm,'BNLRCG') && ~binaryRep)
               error('To use BNLRCG Algorithm, spikeTrain must have a binary representation. Increase sampleRate and try again');
            end
           
           %If performing batchMode analysis, this stacks up the
           %corresponding spike vectors and the design matrices
           for i=1:length(indices)
               if(i==1)
                    y=tObj.getSpikeVector(indices(i));
                    X=tObj.getDesignMatrix(indices(i));
                    lambdaTime = tObj.getCov(1).time;
               else
                    y=[y; tObj.getSpikeVector(indices(i))];
                    X=[X; tObj.getDesignMatrix(indices(i))];
                    offset = max(lambdaTime)+1/tObj.sampleRate;
                    lambdaTime = [lambdaTime; (tObj.getCov(1).time +offset)];
                    
               end
           end
           %For a single neuron given covariates,perform the GLM fit. 
%             
%             if(binaryRep)
%                 distribution = 'binomial';
%                 linkfunction = 'logit';
%             else
%                 distribution = 'poisson';
%                 linkfunction = 'log';
%             end
%             size(X)
%             size(y)
                if(strcmp(Algorithm,'GLM'))
                    distribution = 'poisson';
                    linkfunction = 'log';
                    [b,dev,stats] = glmfit(X,y,distribution, 'link', linkfunction,'constant','off');
                elseif(strcmp(Algorithm,'BNLRCG'))
                    distribution = 'binomial';
                    linkfunction = 'logit';
                    rrflag=0; %ML estimation
                    [b,dev,stats] = bnlrCG(X,y,rrflag);          
                    
                else
                    error('Algorithm not supported!');
                end
                b=real(b); %make sure to avoid complex coefficients ... sometimes algorithms return
                           %complex b with the imaginary part near zero.
                           %Need to explore why. For now just keep the real
                           %part.
                if(length(b)>=1)
                    if(strcmp(distribution,'binomial'))
                        data = exp(X*b(1:end));
                        data = (data./(1+data)).*tObj.sampleRate;
                        
                    elseif(strcmp(distribution,'poisson'));
                        data = exp(X*b(1:end)).*tObj.sampleRate;
                        
%                         
                    end
                end
                

                
                lambdaIndexStr = num2str(lambdaIndex); 
                
                lambda=Covariate(lambdaTime,data,...
                       '\lambda(t)',tObj.getCov(1).xlabelval,...
                        tObj.getCov(1).xunits,'Hz',strcat('\lambda_{',lambdaIndexStr,'}'));
                                mu=b;
                        s=stats.se;
%             Mc=30;
%             for c=1:Mc
%                 z=normrnd(0,1,length(s),1);
%                 bKDraw(:,c)=mu+(s.*z);
%             end
%             if(strcmp(distribution,'poisson'))
%                 lambdaDraw=exp(X*bKDraw)*(tObj.sampleRate);
%             else 
%                 lambdaDraw=exp(X*bKDraw)./(1+exp(X*bKDraw))*(tObj.sampleRate);
%             end
%             lambdaDraw(isinf(lambdaDraw))=0;
%             alphaVal=.05;
%             for k=1:length(lambdaDraw)
%               [f,x] = ecdf(squeeze(lambdaDraw(k,:)));
%               CIs(k,1) = x(find(f<alphaVal/2,1,'last'));
%               CIs(k,2) = x(find(f>(1-alphaVal/2),1,'first'));
%             end
%            
%   
%             ciPSTHGLM = ConfidenceInterval(lambdaTime,CIs,'CI_{psth_GLM}',lambda.xlabelval,lambda.xunits,lambda.yunits);
%             lambda.setConfInterval(ciPSTHGLM);
            
                %The deviance should be real since it a probability measure
                %and therefore any imaginary part is ignored.
                AIC = 2*length(b)+real(dev);
                BIC = length(b)*log(length(y))+real(dev);
        end        
        function handle = plotInvGausTrans(fitResults,makePlot)
            % handle = plotInvGausTrans(fitResults,makePlot)
            % Given the CDF of the rescaled spike times (the u'js) computes
            % the auto-correlation function inverse gaussian tranformated
            % u'js and the 95% confidence interval that they are distinct
            % from zero. 
            % Idea: if gaussian RVs are uncorrelated, they are indep., then
            %       this suggest independence of the uj's and of the zj's
            %       from the time-rescaling theorem. If zj's are
            %       independent and KS plot is within 95% confidence
            %       interval suggests that candidate lambda is close to the
            %       true lambda.
            if(nargin<2)
                makePlot=0;
            end
            [X,rhoSig,confBoundSig] = Analysis.computeInvGausTrans(fitResults.Z);
            fitResults.setInvGausStats(X,rhoSig,confBoundSig);
            
            if(fitResults.isValDataPresent)
                  [X,rhoSig,confBoundSig] = Analysis.computeInvGausTrans(fitResults.validation.Z);
                  fitResults.validation.setInvGausStats(X,rhoSig,confBoundSig);
            end
            
            if(makePlot==1)
                handle=fitResults.plotInvGausTrans;
            end
            
        end
        function plotFitResidual(fitResults,windowSize,makePlot)
            % plotFitResidual(fitResults,windowSize,makePlot)
            % computes the point process residual between the true spike
            % train and that predicted by the candidate conditional
            % intensity function.
            % The result is stored in fitResult.
            %
            if(nargin<3 || isempty(makePlot))
                makePlot=1;
            end
            if(nargin<2 || isempty(windowSize))
                windowSize=.01;
            end
            M = Analysis.computeFitResidual(fitResults.neuralSpikeTrain,fitResults.lambda,windowSize);
            fitResults.setFitResidual(M);
            
            if(fitResults.isValDataPresent)
                 M = Analysis.computeFitResidual(fitResults.validation.neuralSpikeTrain,fitResults.validation.lambda,windowSize);
                 fitResults.validation.setFitResidual(M);
            end
 
            if(makePlot)
              fitResults.plotResidual;
            end
        end
        function handle = KSPlot(fitResults,DTCorrection,makePlot)
            %handle = KSPlot(fitResults,makePlot)
            % Computes the KS statistics and makes the plot. Stores
            % appropriate parameters in fitResults.
            % If validation data is also available, it does the same for
            % the validation data.
            % DTCorrection: 0 for no DT Correction of KS plot, 1 is the
            %               default.
            if(nargin <3)
                makePlot =1; %By default make the plot
            end
            if(nargin<2)
               DTCorrection = 1; 
            end
                

            [Z, U, xAxis, KSSorted, ks_stat] = Analysis.computeKSStats(fitResults.neuralSpikeTrain,fitResults.lambda,DTCorrection);
            fitResults.setKSStats(Z,U, xAxis, KSSorted, ks_stat);
            
            
            if(fitResults.isValDataPresent)
                 %make sure nst is in appropriate window
                    [Z, U, xAxis, KSSorted, ks_stat] = Analysis.computeKSStats(fitResults.validation.neuralSpikeTrain,fitResults.validation.lambda,DTCorrection);
                    fitResults.validation.setKSStats(Z, U, xAxis, KSSorted, ks_stat);
   
            end
            
            if(makePlot)
                handle = fitResults.KSPlot; 
            else
                handle = [];
            end
                
        end  
        function handle = plotSeqCorr(fitResults)
            % handle = plotSeqCorr(fitResults)
            % Plots the sequential correlation coefficients of the rescaled
            % ISIs. zj vs. zj-1
            handle = fitResults.plotSeqCorr;
            
        end
        function handle = plotCoeffs(fitResults)
            % handle = plotCoeffs(fitResults)
            % Plots the regression coefficients for all the different fits.
            
            handle = fitResults.plotCoeffs;
            
        end
        
        
        function [X,rhoSig,confBoundSig]   = computeInvGausTrans(Z)
            % [U,X,rhoSig,confBoundSig]   = computeInvGausTrans(Z)
            % Take rescaled spikeTimes, zjs, transforms them to
            % uniform(0,1), then computes the inverse gaussian 
            % transformation of these to xj's. rhoSig is the
            % auto-correlation funcion of these xj's and is used to test
            % for independence of the xj's. Independence of the xj's
            % suggests indepence of the uj's and zj's (a condition
            % necessary for the Time Rescaling Theorem).
            
            U=1-exp(-Z);
            U(U>=.999999)=.999999; %Prevent any 1 values which lead to infinity in X
            U(U==0)=.000001;
            X = norminv(U,0,1);
            %X=erfinv(U);
            [~,colm] = size(X);
            for i=1:colm
                [c(:,i),lags] = xcov(X(:,i));
            end
            index=find(lags==1);
            lags=lags(index:end);
            rho=c(index:end,:)./repmat(c(index-1,:),length(lags),1);
            n=length(X);
            % Defaults to the 95% confidence intervals
            % Can extend to allow selection of 95% or 99% CI
            confBound = 1.96/sqrt(n)*ones(length(lags),1);
%              size(lags)
%              size(rho)
            
            confBoundSig = SignalObj(lags,[confBound -confBound],'ACF[ \Phi^{-1}(u_i) ]','\Delta \tau','sec');
            confBoundSig.setPlotProps({' ''r'', ''LineWidth'' ,3'},1);
            confBoundSig.setPlotProps({' ''r'', ''LineWidth'' ,3'},2);

            handle=[];
            rhoSig = SignalObj(lags,rho, 'ACF[ \Phi^-1(u_i) ]','Lag \Delta \tau','sec');
            plotProps = cell(1,colm);
            for i=1:colm
                plotProps{i}=strcat('''', '.',Analysis.colors{mod(i-1,length(Analysis.colors))+1},'''');
            end
            rhoSig.setPlotProps(plotProps);
   
            
            
        end        
        function [Z,U,xAxis,KSSorted, ks_stat] = computeKSStats(nspikeObj,lambdaInput,DTCorrection)
            % [Z,U,xAxis,KSSorted, ks_stat] = computeKSStats(nspikeTrain,lambdaInput)
            % Given a neural spike train (a sequence of spike times) and a
            % conditional intensity function, computes the rescaled ISIs
            % according to the time-rescaling theorem in Z. The Uj are
            % returned in U and correspond to a transformation fo the Zjs
            % (exponential rate 1 (according to T-R theorem) to be
            % uniform(0,1). 
            %
            % DTCorrection: 0 for no DT Correction of KS plot, 1 is the
            %               default.
            % nspikeTrain: a nspikeTrain object
            % lambdaInput: candidate conditional intensity function (a Covariate)
            % Z: rescaled spike times
            % U: Zjs tranformed to be uniform(0,1) 
            % xAxis: x-axis of the KS plot
            % KSSorted: y-axis of KS plot
            % ks_stat: the KS statistic. Maximum deviations from the 45
            % degree line for each conditional intensity function.
            
            %get the relevant spike train
            if(nargin<3)
                DTCorrection =1; 
            end
            

            if(length(nspikeObj)>1) %in batch analysis we get multiple trials
                nstCollObj = nstColl(nspikeObj);
                nCopy = nstCollObj.toSpikeTrain;
               
            else
                nCopy =nspikeObj.nstCopy;
                
            end
            
  
%             minTime = nCopy.minTime;
%             maxTime = nCopy.maxTime;
            nCopy.resample(lambdaInput.sampleRate);
            nCopy.setMinTime(lambdaInput.minTime);
            nCopy.setMaxTime(lambdaInput.maxTime);
            
            repBin = nCopy.isSigRepBin; 
            if(~repBin)
               lambdaInput=lambdaInput.resample(2*lambdaInput.sampleRate);
               nCopy.resample(lambdaInput.sampleRate);
            end

            if(DTCorrection==1 && repBin)  
                % Use DT Correction for Time Rescaling Theorem - Haslinger, Pipa and Brown (2010)
                pkSignal=lambdaInput;
                pk = pkSignal.data.*(1/lambdaInput.sampleRate);
                pk = max(pk,1e-10);
                spikeTrain = nCopy.getSigRep.data;
                minDim = min(size(pk,1),size(spikeTrain,1));
                pk=pk(1:minDim,:);
                spikeTrain=spikeTrain(1:minDim,:);
                
                intValues=zeros(length(nCopy.getSpikeTimes)-1,lambdaInput.dimension);
                for i=1:lambdaInput.dimension
                    pk(:,i) = nanmin(nanmax(pk(:,i),0),1);
                    temp = ksdiscrete(pk(:,i),spikeTrain,'spiketrain');
%                     length(temp)
%                     length(intValues(:,i))
                    %sometimes ksdiscrete returns 1 less spike train than
                    %expected ... need to debug .... for now just fix
                    %using length(temp) to index into intValues;
                    intValues(1:length(temp),i) = temp;
                end

                
            else % do not correct for discrete time effects
                
                tempLambda = lambdaInput;
%                 tempLambda = tempLambda.resample(tempLambda.sampleRate*4);
%                 lambda=tempLambda.getSigInTimeWindow(minTime,maxTime);%.dataToMatrix;
                lambdaPosdata = max(tempLambda.data,0);
                lambda = Covariate(tempLambda.time,lambdaPosdata,tempLambda.name,tempLambda.xlabelval,tempLambda.xunits,tempLambda.yunits,tempLambda.dataLabels);
                lambdaInt = lambda.integral;
                

                if(nCopy.isSigRepBin)
                    spikeTimes = nCopy.getSpikeTimes;
                    spikeTimes = [0 spikeTimes];
                    
                else
%                     spikeTimes = nCopy.getSpikeTimes;
%                     maxBinSize=nCopy.getMaxBinSizeBinary;
%                     lambdaInt = lambda.resample(1/maxBinSize).integral;
                     nstSignal = nCopy.getSigRep;
                     spikeTimes=nstSignal.time(nstSignal.data~=0);
                     spikeTimes = [0 spikeTimes'];

                end

                   if(~isempty(spikeTimes))
                        tempVals = lambdaInt.getValueAt(spikeTimes);                        
                        intValues= tempVals(2:end,:)-tempVals(1:end-1,:);
                    else
                        intValues = 0;
                   end

    %                 intValues=2*intValues;
    %             lambdaInt.plot; hold all;
    %             vals =lambdaInt.getValueAt(spikeTimes);
    %             plot(spikeTimes,vals,'.')
                
            end
            Z = intValues; % rescales spike times - exponential rate 1
            U = 1-exp(-Z); % store the rescaled spike times - uniform(0,1)
             

            KSSorted = sort( U,'ascend' );
            N = size(KSSorted,1);
            if(N~=0)
                xAxis=(([1:N]-.5)/N)'*ones(1,lambdaInput.dimension);
                ks_stat = max(abs(KSSorted - (([1:N]-.5)/N)'*ones(1,lambdaInput.dimension))); 
            else
                ks_stat=1;
                xAxis=[];
            end
        end
        function M=computeFitResidual(nspikeObj,lambda,windowSize)
            % M=computeFitResidual(nspikeTrain,lambda,windowSize)
            % Computes the Point Process residual defined in
            % 'A point process framework for relating neural spiking
            % activity to spiking history, W Truccolo, UT Eden, MR Fellows,
            % JP Donoghue and EN. Brown. Journal of Neurophysiology 2005.
            %
            % nspikeTrain: nspikeTrain object
            % lambda: candidate conditional intensity function evaluated on the time
            %         interval of the spike train.
            % windowSize: the size of the window over which to compute the
            %             residual.
            % M: the point process residual (a Covariate object).
            %
            if(nargin<3 || isempty(windowSize))
                windowSize=.01;
            end
            

            if(length(nspikeObj)>1) %in batch analysis we get multiple trials
                nstCollObj = nstColl(nspikeObj);
                nCopy = nstCollObj.toSpikeTrain;
            else
                nCopy =nspikeObj.nstCopy;
            end
            
            nCopy.resample(lambda.sampleRate);
            nCopy.setMinTime(lambda.minTime);
            nCopy.setMaxTime(lambda.maxTime);
            
            
            sumSpikes=nCopy.getSigRep(windowSize);%tObj.getNeuron(fitResults.neuronNumber).nstCopy;
%             sumSpikesOverWindow = sumSpikes.data(1:end);
%             windowTimes = nCopy.minTime:windowSize:lambda.maxTime;
            windowTimes = linspace(nCopy.minTime,nCopy.maxTime,length(sumSpikes.time));           
            lambdaInt = lambda.integral;
            lambdaIntVals = lambdaInt.getValueAt(windowTimes(2:end))-lambdaInt.getValueAt(windowTimes(1:(end-1)));
            if(length(lambdaIntVals)==length(sumSpikes.data))
                sumSpikesOverWindow = sumSpikes.data(1:end);
            elseif(length(lambdaIntVals)<length(sumSpikes.data))
                sumSpikesOverWindow = sumSpikes.data(2:end);
            end
            Mdata=repmat(sumSpikesOverWindow,[1 lambdaInt.dimension])-lambdaIntVals;
            dataLabels = cell(1,lambdaInt.dimension);
            for i=1:lambdaInt.dimension
                dataLabels{i} = lambda.dataLabels{i};
            end
            
            M=Covariate(windowTimes(1:end),[zeros(1,size(Mdata,2));Mdata],'M(t_k)',lambdaInt.xlabelval, ...
                        lambdaInt.xunits,lambdaInt.yunits,dataLabels);
            
        end
        
        function [fitResults,ensembleCovariate,tcc] = compHistEnsCoeffForAll(tObj,history,makePlot)
            % [fitResults,ensembleCovariate,tcc] = compHistEnsCoeffForAll(tObj,history,makePlot)
            %  runs Analysis.compHistEnsCoff for each neuron that is not masked.
            if(nargin<3 || isempty(makePlot))
                makePlot=1;
            end
            neuronIndex=tObj.getNeuronIndFromMask;
            fitResults = cell(1,length(neuronIndex));
            tcc = cell(1,length(neuronIndex));
            ensembleCovariate = tObj.getEnsembleNeuronCovariates(neuronIndex(1),[],history);
            [fitResults{1},tcc{1}] = Analysis.compHistEnsCoeff(tObj,history,neuronIndex(1),tObj.getNeuronNeighbors(1),ensembleCovariate,makePlot);
            for i=2:length(neuronIndex)
                ensembleCovariate.maskAwayAllExcept(tObj.getNeuronNeighbors(neuronIndex(i)));
                [fitResults{i},tcc{i}] = Analysis.compHistEnsCoeff(tObj,history,neuronIndex(i),tObj.getNeuronNeighbors(neuronIndex(i)),ensembleCovariate,makePlot);           
            end
        end
        function [fitResults,ensembleCov,tcc] = compHistEnsCoeff(tObj,history,neuronNum,neighbors,ensembleCov,makePlot)
            % [fitResults,ensembleCov,tcc] = compHistEnsCoeff(tObj,history,neuronNum,neighbors,ensembleCov,makePlot)
            % Given a trial, a history object compute the history time
            % series for the ensemble of neighboring neurons. This is done for all neurons and the result is returned in 
            % ensembleCov as a covariate collection. This collection is
            % then used as the design matrix and the analysis performed for
            % each neuron. The results are returned in fitResults.
            %
            %
            if(nargin<6 || isempty(makePlot))
                makePlot=1;
            end
            
            if(nargin<3 || isempty(neuronNum))
                neuronNum=tObj.getNeuronIndFromMask;
            end
            
            if(nargin<4 || isempty(neighbors))
                 neighbors=tObj.getNeuronNeighbors(neuronNum); %every other neuron
            end
            
            if(nargin<5 || isempty(ensembleCov))
                ensembleCov = tObj.getEnsembleNeuronCovariates(neuronNum,neighbors,history);
            end
            

            %Create a covariate collection that consists of only the
            %ensemble history
            ensembTrial = Trial(tObj.nspikeColl,ensembleCov);
            tc=TrialConfig('all',[],[]); %use all ensembleCov
            tcc = ConfigColl(tc); 
            fitResults =Analysis.RunAnalysisForNeuron(ensembTrial,neuronNum,tcc,makePlot);
        end
        function [fitResults,tcc] = computeHistLag(tObj,neuronNum,windowTimes,CovLabels,Algorithm,batchMode,sampleRate,makePlot,histMinTimes,histMaxTimes)
            % [fitResults,tcc] = computeHistLag(tObj,tObj,neuronNum,windowTimes,CovLabels,sampleRate,makePlot)
            % For the neuron in neuronNum, runs an analysis with self
            % history, and no extrinsic covariates, and no ensemble history
            % as a covariates. The self history is specfied by a vector
            % of windowTimes. There will be length(windowTimes) different
            % results (configurations) corresponding to increasing number of history
            % windows.
            if(nargin<10)
                histMaxTimes =[];
            end
            if(nargin<9)
                histMinTimes=[];
            end
            if(nargin<8)
                makePlot=1;
            end
            if(nargin<7 || isempty(sampleRate))
                sampleRate = tObj.sampleRate;
            end
            if(nargin<6 || isempty(batchMode))
                batchMode = 0;
            end
            if(nargin<5 || isempty(Algorithm))
                Algorithm = 'GLM';
            end
            if(nargin<4)
                CovLabels ={};
            end
            if(nargin<3)
                error('Must specify a vector of windowTimes');
            end
            if(nargin<2 || isempty(neuronNum))
                neuronNum = tObj.getNeuronIndFromMask;
            end
            
            % tcObj=TrialConfig(covMask,sampleRate, history,minTime,maxTime)
            tc=cell(1,length(windowTimes)-1);
            for i=1:length(tc)+1
                %use no covariates
                if(i==1)
                    tc{i} = TrialConfig(CovLabels,sampleRate,[],[]); tc{i}.setName('Baseline');
                else
                    if(and(isempty(histMinTimes),isempty(histMaxTimes)))
                        tc{i} = TrialConfig(CovLabels,sampleRate,windowTimes(1:i)); tc{i}.setName(strcat('Window',num2str(i-1)));
                    else
                        hTemp = History(windowTimes(1:i),histMinTimes,histMaxTimes);
                        tc{i} = TrialConfig(CovLabels,sampleRate,hTemp); tc{i}.setName(strcat('Window',num2str(i-1)));
                    end
                end
                    
            end
            DTCorrection=1;
            tcc = ConfigColl(tc);             
            
            fitResults =Analysis.RunAnalysisForNeuron(tObj,neuronNum,tcc,makePlot,Algorithm,DTCorrection,batchMode);
                                 
        end
        function fitResults = computeHistLagForAll(tObj,windowTimes,CovLabels,Algorithm,batchMode,sampleRate,makePlot,histMinTimes,histMaxTimes)
            % [fitResults,tcc] = computeHistLagAll(tObj,windowTimes,CovLabels,sampleRate,makePlot)
            % Calls computeHistLab for each neuron in the trial that is not masked. 
            if(nargin<9)
                histMaxTimes =[];
            end
            if(nargin<8)
                histMinTimes=[];
            end
            if(nargin<7)
                makePlot=1;
            end
            if(nargin<6 || isempty(sampleRate))
                sampleRate = tObj.sampleRate;
            end
            if(nargin<5 || isempty(batchMode))
                batchMode=0;
            end
            if(nargin<4 || isempty(Algorithm))
                Algorithm = 'GLM';
            end
            
            if(nargin<3)
                CovLabels ={};
            end
            if(nargin<2)
                error('Must specify a vector of windowTimes');
            end
            
            neuronIndex=tObj.getNeuronIndFromMask;
            fitResults = cell(1,length(neuronIndex));
            for i=1:length(neuronIndex)
               fitResults{i} = Analysis.computeHistLag(tObj,neuronIndex(i),windowTimes,CovLabels,Algorithm,batchMode,sampleRate,makePlot,histMinTimes,histMaxTimes);
            end
            

            
        end
        function [fitResults,tcc]=computeNeighbors(tObj,neuronNum,sampleRate,windowTimes,makePlot)
            % [fitResults,tcc]=computeNeighbors(tObj,neuronNum,sampleRate,windowTimes,makePlot)
            % For the neuron in neuronNum, runs an analysis with no self
            % history, and no extrinsic covariates, only ensemble history
            % as a covariate. The ensemble history is specfied by a vector
            % of windowTimes. There will be length(windowTimes) different
            % results corresponding to increasing number of history
            % windows.
            if(nargin<4)
                makePlot=1;
            end
            if(nargin<3 || isempty(sampleRate))
                sampleRate = tObj.sampleRate;
            end
            if(nargin<2 || isempty(neuronNum))
                neuronNum = tObj.getNeuronIndFromMask;
            end
            tc=cell(1,length(windowTimes)-1);
            for i=1:length(windowTimes)
                % For reference: tcObj=TrialConfig(covMask,sampleRate, history,covHist,covLag,name)
                if(i==1)
                    tc{i} = TrialConfig({},sampleRate,[],[]); tc{1}.setName('Baseline');
                else
                    tc{i} = TrialConfig({},sampleRate,[],windowTimes(1:i)); 
                end
            end
               tcc = ConfigColl(tc);             
               fitResults =Analysis.RunAnalysisForNeuron(tObj,neuronNum,tcc,makePlot);
        end
        
        function cc=spikeTrigAvg(tObj,neuronNum,windowSize)
            % cc=spikeTrigAvg(tObj,neuronNum,windowSize)
            % Each covariate dimension is sampled at every spike time of
            % the neuron specified +/- windowSize. The number of columns of
            % each new covariate corresponds to the number of spikes in the
            % spike train. A covariate collection is returned containing
            % each covariate dimension as a new covariate. These can then
            % be easily avergaged by using SignalObj method
            % plotVariability.
           
            nCopy=tObj.getNeuron(neuronNum).nstCopy;
            t=-windowSize/2:1/tObj.sampleRate:windowSize/2;
            covIndex=0;
            for i = 1:tObj.covarColl.numCov
                tempCov=tObj.getCov(i);
                data=[];
                dataIndex=0;
 %               for n=1:tObj.nspikeColl.numSpikeTrains
%                     nCopy=tObj.getNeuron(neuronNum).nstCopy;
                    spikeTimes = nCopy.getSpikeTimes';
                    for j=1:length(spikeTimes)
                      dataIndex=dataIndex+1;
                      tc=tempCov.getSigInTimeWindow(spikeTimes(j)-windowSize/2,spikeTimes(j)+windowSize/2);
                      tc=tc.shift(-tc.minTime-windowSize/2);
                      tempData = tc.getValueAt(t);
%                       if(isempty(data))
%                         data(dataIndex,1:length(tempData),:)=tempData;
%                       else
%                         data(dataIndex,:,:)=zeros(size(squeeze(data(1,:,:))));
                        data(dataIndex,:,:)=tempData;
%                       end
                    end
                %end
                
                for k=1:tempCov.dimension
                    covIndex = covIndex+1;
                    cov{covIndex} = Covariate(t,squeeze(data(:,:,k)),tempCov.dataLabels{k},tempCov.xlabelval,tempCov.xunits,tempCov.yunits,tempCov.dataLabels{k}); 
                end
            end
            cc=CovColl(cov);
            
        end
 

    end
    
    
    
end



function [flatMask, maxIndex]=flatMaskCellToMat(flatMaskCell)
    % [flatMask, maxIndex]=flatMaskCellToMat(flatMaskCell)
    lMask =zeros(1,length(flatMaskCell));
    for i=1:length(flatMaskCell)
        lMask(i) = length(flatMaskCell{i});
    end
    [maxSize,maxIndex] = max(lMask);
    flatMask = zeros(maxSize,length(flatMaskCell));
    for i=1:length(flatMaskCell)
        flatMask(1:length(flatMaskCell{i}),i) = flatMaskCell{i};
    end
end
function [beta_new,devnew,stats] = bnlrCG(X,yframe,rrflag)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                               %
% Truncated, L-2 Regularized, Binomial Logistic Regression with % 
% Conjugate Gradient Solver                                     %
%                                                               %
% Author: Demba Elimane Ba                                      %
%         MIT Department of EECS                                %
%         Neuro Stat Research Lab (MIT Department of BCS)       %
% Date  : August the 25th, 2008                                 %
%                                                               %
% Note  : Matlab implementation of Paul Komarek's TR-IRLS       %
%                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Modified by Iahn Cajigas 9-23-09 to automatically add the DC term for the
% design matrix

% Arguments:
%   X:      design matrix, including DC column of all ones (1st or last)
%   yframe: column vector of binary observations
%   rrflag: 1: MAP estimation with Gaussian(0,sigma^2) prior (default value 
%               1/sigma^2 = 10) 
%           0: ML estimation (1/sigma^2 = 0, i.e. sigma -> infinity)
%
%   N.B: The equivalent call with glmfit is:
%
%   [beta,dev,stats]=glmfit(X(:,2:end),[yframe ones(size(yframe))],'binomial','logit');
%
%   Unlike glmfit, this function assumes that the design matrix
%   already has a column of all ones (1st or last).

    %Modify the design Matrix;
    [rows,~]=size(X);
%     X = [ones(rows,1), X];
    % CG parameters
    cgmax = 30;
    cgeps = 1e-6;

    % LR parameters
    lrmax = 100;
    lreps = 0.05;
    lambda = 10;

    [n,d] = size(X);
    % Perform logistic regression
    i = 0;
    % Initial guess for Beta = beta_old ?
    beta_old = zeros(d,1);
    n = X*beta_old;
    u = exp(n)./(1+exp(n));
    W = repmat(u'.*(1-u)',d,1);
    z = X*beta_old + (W(1,:)'.^-1).*(yframe - u);

    devold = -2*sum(yframe.*log(u) + (1-yframe).*log(1-u));
    devnew = 0;
    devdiff = abs(devnew - devold);

% B=[];
    while (i < lrmax && devdiff > lreps)
        % Do CG -> beta_new, i.e. solve for beta_new: XtWX*beta_new = XtWz(beta_old) using CG

        A = X'.*W*X; b = X'.*W*z;

        % A needs to be positive definite so any zero or negative 
        % eigenvalues are set to machine precision.
        [vec,val]=eig(A); val(val<=0)=eps;
        A=vec*val*vec';
                    
        if(any(isnan(b)))
           b(isnan(b))=0; 
        end
        if rrflag == 1
            A = A + lambda*eye(size(A));
        end

        [beta_new, flag] = cgs(A,b,cgeps,cgmax,[],[],beta_old);    
        beta_old = beta_new;

        n = X*beta_old;
        u = exp(n)./(1+exp(n));
        W = repmat(u'.*(1-u)',d,1);
        z = X*beta_old + (W(1,:)'.^-1).*(yframe - u);

        devnew = -2*sum(yframe.*log(u) + (1-yframe).*log(1-u));
        devdiff = abs(devnew - devold);    
        devold = devnew;

        i = i+1;

%         B=[B,beta_new];
    end

    % Compute a few statistics
    stats.dfe = 0;
    stats.s = 0;
    stats.sfit = 0;
    stats.covb = inv(A);
    stats.se = sqrt(diag(stats.covb));
    stats.coeffcorr = stats.covb./sqrt((repmat(diag(stats.covb),1,d).*repmat(diag(stats.covb)',d,1)));
    stats.t = 0;
    stats.p = 0;
    stats.resid = 0;
    stats.residp = 0;
    stats.residd = 0;
    stats.resida = 0;
end


function [rst,varargout] = ksdiscrete(pk,st,spikeflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ksdiscrete.m 
% written by Rob Haslinger, December 2009
%
% This function performs time rescaling of ISIs based upon the discrete
% time version of the time rescaling theorem as described in Haslinger,
% Pipa and Brown (2010?).  This method corrects for biases in the KS plot 
% caused by the temporal discretization.
%
% This function can be called in two ways
%
% 1) input the discrete time conditional probabilities "pk"  where 0<=pk<= 1
% and the spike train "spiketrain" which has elements either equal to 0 (no
% spike) or 1 (spike). There is also a flag 'spiketrain' to indicate that
% it is the full spike train.
%
% [rst,rstsort,xks,cb,rstoldsort] = ksdiscrete(pk,spiketrain,'spiketrain')
%
% 2) input the discrete time conditional probabilities "pk" and a list of
% the indicies "spikeind" of the bin indicies that the spikes are locaed in. 
% There is also a flag 'spikeind' to indicate that the indicies are
% being given, not the full spike train
%
% [rst,rstsort,xks,cb,rstoldsort] = ksdiscrete(pk,spikeind,'spikeind');
%
% required output:
%
% rst : a vector of unsorted uniformly distributed rescaled times. This is
% the only output that is required.
% 
% optional output, given in the order they appear in the list function
% outputs :
%
% rstsort : a vector of rescaled times sorted into ascending order
% xks : a vector of x axis values to plot the sorted rescaled times against
% cb : the value of the 95% confidence bounds
% rstoldsort : a vector of sorted rescaled times done without the discrete
% time correction
%
% To make a KS plot one would do
% plot(xks,rstsort,'k-');
% hold on;
% plot(xks,xks+cb,'k--',xks,xks-cb,'k--');
%
% To make a Differential KS plot one would do
% plot(xks,rstsort-xks,'k-');
% hold on;
% plot(xks,zeros(length(xks))+cb,'k--',xks,zeros(length(xks))-cb);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with determining the inputs and some basic input error checking

    if nargin < 3 || nargin > 3;  
        error('Number of input arguments must be equal to 3'); 
    end;

    % make pk into a column vector;

    [m1,m2]=size(pk);
        if (m1 ~=1 && m2 ~=1); error('pk must be a vector'); end;
        if (m2>m1); pk=pk'; end;
        [m1,m2]=size(pk);

    % make sure pk's are within [0,1]
    index=find(pk<0);
    if isempty(index) ~=1; 
        error('all values for pk must be within [0,1]'); 
    end;
    index=find(pk>1);
    if isempty(index) ~=1; 
        error('all values for pk must be within [0,1]'); 
    end;
    clear index;    

    % make column vector of spike indicies

    if strcmp(spikeflag,'spiketrain'); % spike train input

        [n1,n2]=size(st);
          if (n1 ~=1 && n2 ~=1); error('spike train must be a vector'); end;
        if (n2>n1); st=st'; end;

        if m1 ~= n1; error('pk and spike train must be same length'); end;

        spikeindicies=find(st==1);

        Nspikes=length(spikeindicies);

    elseif strcmp(spikeflag,'spikeind'); % spike index input

        [n1,n2]=size(st);
          if (n1 ~=1 && n2 ~=1); error('spike indicies must be a vector'); end;
        if (n2>n1); st=st'; end;

        spikeindicies=unique(st);
        Nspikes=length(spikeindicies);

    end;

    % check that those indicies are in [1:length(pk)];

    if(isempty(spikeindicies))
    	rst = pk;
        return;
    end
    if spikeindicies(1)<1; 
         error('There is at least one spike with index less than 0'); 
    end;
    if spikeindicies(Nspikes)>length(pk); 
         error('There is at least one spike with a index greater than the length of pk'); 
    end;    

    % error checking done

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now do the actual discrete time KS test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % initialize random number generator
    s = RandStream('mt19937ar','Seed', sum(100*clock));
    RandStream.setGlobalStream(s);
    %rand('twister',sum(100*clock));

    % make the qk's

    qk=-log(1-pk);

    % make the rescaled times

    rst=zeros(Nspikes-1,1);
    rstold=zeros(Nspikes-1,1);

    for r=1:Nspikes-1;

        total = 0;

        ind1=spikeindicies(r);
        ind2=spikeindicies(r+1);

        total=total+sum(qk(ind1+1:ind2-1));

        delta=-(1/qk(ind2))*log(1-rand()*(1-exp(-qk(ind2))));

        total=total+qk(ind2)*delta;

        rst(r)=total;

        rstold(r)=sum(qk(ind1+1:ind2));

    end;

%     rst=1-exp(-rst);
%     rstold=1-exp(-rstold);

    % optional outputs

    rstsort=sort(rst);
    varargout{1}=rstsort;

    inrst=1/(Nspikes-1);
    xrst=(0.5*inrst:inrst:1-0.5*inrst)';
    varargout{2}=xrst;

    cb=1.36*sqrt(inrst);
    varargout{3}=cb;    

    varargout{4}=sort(rstold);
end
