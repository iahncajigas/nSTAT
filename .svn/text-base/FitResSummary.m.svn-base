classdef FitResSummary < handle
% FITRESSUMMARY - given a collection of fitResult objects (one for each
% neuron, each containing the results of multiple regressions), computes
% summary statistics across neurons. This is to allows visualization of
% commonalities in the data across multiple neurons. 
%
% % <a href="matlab: methods('FitResSummary')">methods</a>
% see also <a href="matlab:help('FitResult')">FitResult</a>
%
% Reference page in Help browser
% <a href="matlab:doc('FitResSummary')">doc FitResSummary</a>    
    

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
        fitResCell %collection of FitResult Objects
        fitNames   %names of the different fits
        numResults %number of fits in each FitResult Object
        maxNumIndex %Index into the results that have all of the labels
        dev
        AIC % AIC matrix (length(neuronNumbers) x numResults)
        BIC % BIC matrix (length(neuronNumbers) x numResults)
        bAct
        seAct
        sigIndex
        covLabels %labels of relevant covariates in fits.
        uniqueCovLabels
        indicesToUniqueLabels
        flatMask
        neuronNumbers % neuronNumber(i) is the number of the neuron corresponding to fitResCell{i}
        numNeurons  % Total number of neurons being summarized
        numCoeffs 
        numResultsCoeffPresent;
        KSStats     % KSStats matrix (length(neuronNumbers) x numResults)
        KSPvalues
        withinConfInt
        plotParams
        coeffRange  
     end
    
    methods 
        function frsObj = FitResSummary(fitResultsCell)
            % frsObj = FitResSummary(fitResultsCell)
            % Takes a cell array of FitResult objects and returns a
            % FitResultSummary object.
           if(isa(fitResultsCell,'FitResult'))
               frsObj=FitResSummary({fitResultsCell});
           elseif(isa(fitResultsCell,'cell') && ~isempty(fitResultsCell))
               if(isa(fitResultsCell{1},'FitResult'))
                   frsObj.fitResCell=cell(length(fitResultsCell),1);
                   maxNumResults = 0;
                   maxNumIndex=1;
                   for i=1:length(fitResultsCell)
                      if(fitResultsCell{i}.numResults>maxNumResults)
                          maxNumResults = fitResultsCell{i}.numResults;
                          maxNumIndex   = i;
                      end

                   end    
                   frsObj.maxNumIndex   = maxNumIndex;
                   frsObj.dev           = nan(length(fitResultsCell),maxNumResults);
                   frsObj.AIC           = nan(length(fitResultsCell),maxNumResults);
                   frsObj.BIC           = nan(length(fitResultsCell),maxNumResults);
                   frsObj.KSStats       = nan(length(fitResultsCell),maxNumResults);
                   frsObj.KSPvalues     = nan(length(fitResultsCell),maxNumResults);
                   frsObj.withinConfInt = zeros(length(fitResultsCell),maxNumResults);
                   for i=1:length(fitResultsCell)
                        frsObj.neuronNumbers(i)           = fitResultsCell{i}.neuronNumber;
                        frsObj.dev(i,1:length(fitResultsCell{i}.dev))                   = fitResultsCell{i}.dev(1:end);
                        frsObj.AIC(i,1:length(fitResultsCell{i}.AIC))                   = fitResultsCell{i}.AIC(1:end);
                        frsObj.BIC(i,1:length(fitResultsCell{i}.BIC))                   = fitResultsCell{i}.BIC(1:end);
                        frsObj.covLabels{i}               = fitResultsCell{i}.uniqueCovLabels;
                       %flatMask(:,:,i)                   = fitResultsCell{i}.flatMask;
                        frsObj.fitResCell{i}              = fitResultsCell{i};
                        frsObj.KSStats(i,1:length(fitResultsCell{i}.KSStats.ks_stat))               = fitResultsCell{i}.KSStats.ks_stat(1:end);
                        frsObj.KSPvalues(i,1:length(fitResultsCell{i}.KSStats.pValue))               = fitResultsCell{i}.KSStats.pValue(1:end);
                        frsObj.withinConfInt(i,1:length(fitResultsCell{i}.KSStats.withinConfInt))         = fitResultsCell{i}.KSStats.withinConfInt(1:end);
                   end
                   frsObj.numResults = maxNumResults;%fitResultsCell{1}.numResults;
                   frsObj.numNeurons = length(frsObj.neuronNumbers);
                   frsObj.uniqueCovLabels = getUniqueLabels(frsObj.covLabels);
                   frsObj.mapCovLabelsToUniqueLabels;
                   %indActCoeffs = find( sum(squeeze(sum(frsObj.flatMask,2)),2)>0);
                   bAct=nan(length(frsObj.uniqueCovLabels),frsObj.numResults,frsObj.numNeurons);
                   seAct=nan(length(frsObj.uniqueCovLabels),frsObj.numResults,frsObj.numNeurons);
                   sigIndex=zeros(length(frsObj.uniqueCovLabels),frsObj.numResults,frsObj.numNeurons);
                   for n=1:frsObj.numNeurons
                       for j=1:frsObj.numResults
                           index=frsObj.indicesToUniqueLabels{j,n};
                           if(j<=size(fitResultsCell{n}.flatMask,2))
    %                         origIndex = find(fitResultsCell{n}.flatMask(:,j));
                            origIndex = fitResultsCell{n}.indicesToUniqueLabels{j};
                            bAct(index,j,n)     = fitResultsCell{n}.getPlotParams.bAct(origIndex,j);
                            seAct(index,j,n)    = fitResultsCell{n}.getPlotParams.seAct(origIndex,j);
                            sigIndex(index,j,n) = fitResultsCell{n}.getPlotParams.sigIndex(origIndex,j);
                           end
                       end
                   end 
                   frsObj.bAct = bAct;
                   frsObj.seAct= seAct;
                   frsObj.sigIndex=sigIndex;


                   frsObj.numResultsCoeffPresent = (sum(sum(frsObj.flatMask,3),2));               
                   frsObj.numCoeffs  = length(frsObj.uniqueCovLabels);
                   frsObj.coeffRange = [];
                   frsObj.fitNames= fitResultsCell{maxNumIndex}.configNames;

                   frsObj.plotParams.bAct = bAct;%(sum(sum(sigIndex,3),2)>0,:,:);  
                   frsObj.plotParams.seAct = seAct;%(sum(sum(sigIndex,3),2)>0,:,:);  
                   frsObj.plotParams.sigIndex = sigIndex;%(sum(sum(sigIndex,3),2)>0,:,:);  
                   frsObj.plotParams.xLabels  = frsObj.uniqueCovLabels;%(sum(sum(sigIndex,3),2)>0);
                   frsObj.plotParams.numResultsCoeffPresent = frsObj.numResultsCoeffPresent;%(sum(sum(sigIndex,3),2)>0);

    %                frsObj.plotParams.bAct = bAct(sum(sum(sigIndex,3),2)>0,:,:);  
    %                frsObj.plotParams.seAct = seAct(sum(sum(sigIndex,3),2)>0,:,:);  
    %                frsObj.plotParams.sigIndex = sigIndex(sum(sum(sigIndex,3),2)>0,:,:);  
    %                frsObj.plotParams.xLabels  = frsObj.uniqueCovLabels(sum(sum(sigIndex,3),2)>0);
    %                frsObj.plotParams.numResultsCoeffPresent = frsObj.numResultsCoeffPresent(sum(sum(sigIndex,3),2)>0);
    %                
               end
           end
        end

        %% Utility Functions
        function mapCovLabelsToUniqueLabels(frsObj)
            % mapCovLabelsToUniqueLabels(frsObj)
            % from all the covariate labels across all neurons and all
            % fits, finds a minimal list of covariate labels to be used in
            % summarizing the data.
            flatMask = zeros(length(frsObj.uniqueCovLabels),frsObj.numResults,frsObj.numNeurons);
            for n=1:frsObj.numNeurons
                currFitResult =  frsObj.fitResCell{n};
                for j=1:currFitResult.numResults
                    currLabels = currFitResult.covLabels{j};
                    index=zeros(1,length(currLabels));
                        for i=1:length(currLabels)
                            index(i)=strmatch(currLabels{i}, frsObj.uniqueCovLabels, 'exact');
                        end

                        frsObj.indicesToUniqueLabels{j,n} = index;
                        flatMask(index,j,n) = 1;
                end
            end
            frsObj.flatMask = flatMask;
        end
        
        function [dAIC, handle] = getDiffAIC(frsObj,diffIndex,makePlot,h)
            % [dAIC, handle] = getDiffAIC(frsObj,diffIndex,makePlot,h)
            % Takes the AIC matrix and returns a matrix with N-1 columns
            % containing the difference between all columns of the original
            % matrix minus the column indicated by diffIndex. The zero
            % column corresponding to diffIndex is removed from the orginal
            % AIC matrix
            if(nargin<4)
                h=gca;
            end
            if(nargin<3 || isempty(makePlot))
                makePlot=1;
            end
            if(nargin<2 || isempty(diffIndex))
                diffIndex = 1;
            end
            if(frsObj.numResults>1)
                dAIC=computeDiffMat(frsObj.AIC,diffIndex);
            else
                dAIC=frsObj.AIC;
            end
            if(makePlot)
                handle=frsObj.boxPlot(dAIC,diffIndex,h);
            end
        end
        function [dBIC, handle] = getDiffBIC(frsObj,diffIndex,makePlot,h)
            % [dBIC, handle] = getDiffBIC(frsObj,diffIndex,makePlot,h)
            % Takes the BIC matrix and returns a matrix with N-1 columns
            % containing the difference between all columns of the original
            % matrix minus the column indicated by diffIndex. The zero
            % column corresponding to diffIndex is removed from the orginal
            % BIC matrix
            if(nargin<4)
                h=gca;
            end
            if(nargin<3 || isempty(makePlot))
                makePlot=1;
            end
            if(nargin<2 || isempty(diffIndex))
                diffIndex = 1;
            end
            if(frsObj.numResults>1)
                dBIC=computeDiffMat(frsObj.BIC,diffIndex);
            else
                dBIC=frsObj.BIC;
            end
            if(makePlot==1)
                handle=frsObj.boxPlot(dBIC,diffIndex,h);
            end
            
        end     
        function [N,edges,percentSig] = binCoeffs(frsObj,minVal,maxVal,binSize)
            % [N,edges,percentSig] = binCoeffs(frsObj,minVal,maxVal,binSize)
            % Does a histogram of the regression coefficients across all
            % fits. Also returns an indicator of the percentage of times
            % that a coefficient was significantly different than zero out
            % of all the times that it was used in a regression.
            if(nargin<4)
                binSize=.1;
            end
            if(nargin<3)
                if(isempty(frsObj.coeffRange))
                    %v=axis;
                    maxVal=12;%v(4);
                   frsObj.coeffRange.maxVal = maxVal;
                else
                    %if(exists(frsObj.coeffRange.maxVal))
                        maxVal=frsObj.coeffRange.maxVal;
                    %else
                        %maxVal = 12;
                    %end
                end
                
            end
            if(nargin<2)
                if(isempty(frsObj.coeffRange))
                   % v=axis;
                    minVal=-12;%v(3);
                    frsObj.coeffRange.minVal = minVal;
                else
                     %if(exists(frsObj.coeffRange.maxVal))
                        minVal=frsObj.coeffRange.minVal;
                     %else
                     %   minVal=-12;
                     %end
                end
                
            end
            
           edges=(minVal:binSize:maxVal)';
           sigVals = frsObj.plotParams.bAct;%frsObj.plotParams.sigIndex;
           numPlotCoeffs = length(frsObj.plotParams.xLabels);
           numSig = zeros(1,numPlotCoeffs);
           percentSig=zeros(1,numPlotCoeffs);
           sigValArray=[];
           sigGroup=[];
           for i=1:numPlotCoeffs %num coefficients
                tempsigVals = squeeze(sigVals(i,:,:));
                tempsigVals = tempsigVals(squeeze(frsObj.plotParams.sigIndex(i,:,:))==1);
                %sigValArray = [sigValArray;tempsigVals];
                sigGroup    = [sigGroup; repmat(i,[length(tempsigVals),1])];
                Ntemp=histc(tempsigVals,edges);
                numSig(i) = sum(Ntemp);
                [nr,nc] = size(squeeze(sigVals(i,:,:))); %for this coefficient across all fits
                %percentSig(i) = numSig(i)./(nr*nc)*frsObj.numResults./frsObj.plotParams.numResultsCoeffPresent(i);
                percentSig(i) = numSig(i)./frsObj.plotParams.numResultsCoeffPresent(i);
                N(:,i)=Ntemp./numSig(i); %normalize to 1 (pdf)
           end

        end
        function setCoeffRange(frsObj,minVal,maxVal)
            % setCoeffRange(frsObj,minVal,maxVal)
            % Sets the minimum and maximum value for the coeffRange.
            frsObj.coeffRange.minVal=max(-100,minVal);
            frsObj.coeffRange.maxVal=min(100,maxVal);
        end
        function sigValMat = getSigCoeffs(frsObj,fitNum)
            % sigValMat = getSigCoeffs(frsObj,fitNum)
            % if a fitNum is specified, a 2-d matrix (number of rows =
            % number of GLM coefficients for that fit, number of columns =
            % numer of neurons)
            % Otherwise returns a 3-d matrix indicating the significant
            % coefficients across all regression coefficient, fits, and
            % neurons.
            % nUnique params x num regressions x num neurons in size.
            sigValMat = frsObj.plotParams.bAct.*frsObj.plotParams.sigIndex;
            if(nargin==2 && fitNum>0 && fitNum<frsObj.numResults)
                sigValMat = squeeze(sigValMat(:,fitNum,:));
            end
        end
        
        
        %% Plotting Functions
        function handle = plotIC(frsObj,h)
            % handle = plotIC(frsObj,h)
            % Plots the difference in AIC and BIC from baseline (first regression).
            if(nargin<2)
                h(1) = subplot(2,1,1);
                h(2) = subplot(2,1,2);
            end
            
            makePlot=1;
%             [~, h1] = frsObj.getDiffAIC(1,makePlot,h(1)); ylabel('\Delta AIC');
%             [~, h2] = frsObj.getDiffBIC(1,makePlot,h(2)); ylabel('\Delta BIC');
            subplot(2,1,1); h1=frsObj.getDiffAIC(1); ylabel('\Delta AIC'); 
            subplot(2,1,2); h2=frsObj.getDiffBIC(1); ylabel('\Delta BIC'); 
            
            handle = [h1,h2];
        end
        function handle = plotAllCoeffs(frsObj,h,fitNum,plotProps,plotSignificance,subIndex)
            % handle = plotAllCoeffs(frsObj,h)
            % plots the GLM coefficients for each unique covariate across
            % the multiple types of regressions and across neurons.
            
           if(nargin<5 || isempty(plotSignificance))
               plotSignificance = 1;
           end
           
           if(nargin<4 || isempty(plotProps))
               plotProps = [];
           end
           if(nargin<3)
               fitNum = 1:frsObj.numResults;
           end
           if(nargin<2)
               h = gca;
           end
           if(nargin<6 || isempty(subIndex))
              subIndex = 1:size(frsObj.plotParams.bAct,1);
           end
           
           neuronIndex = 1:frsObj.numNeurons;
           bAct = frsObj.plotParams.bAct(subIndex,fitNum,neuronIndex);
           seAct= frsObj.plotParams.seAct(subIndex,fitNum,neuronIndex);
           Xaxis=repmat(1:length(bAct(:,1,1)),[length(bAct(1,:,1)) 1])';
           for i=1:frsObj.numNeurons
                set(gcf,'CurrentAxes',h);
                handle=errorbar(Xaxis,squeeze(bAct(:,:,i)),squeeze(seAct(:,:,i)),'.');%strcat('.',FitResult.colors{mod(i-1,length(FitResult.colors))+1})); 
                hold on;
           end
%             
%            for j=1:length(frsObj.uniqueCovLabels)
%                 for i=1:frsObj.numResults
% %                    set(gcf,'CurrentAxes',h);
% %                   handle=errorbar(Xaxis,squeeze(bAct(:,:,i)),squeeze(seAct(:,:,i)),'.');%strcat('.',FitResult.colors{mod(i-1,length(FitResult.colors))+1})); 
% %                  boxplot(frsObj.KSStats,frsObj.fitNames);
%                     boxplot(squeeze(bAct(i,j,:)))
%                     hold all;
%                 end
%            end
%            
                      
           hy=ylabel('Fit Coefficients','Interpreter','none');
           xtickLabels = frsObj.plotParams.xLabels(subIndex);
           xticks = 1:(length(xtickLabels));
           
           set(gca,'xTick',xticks,'xTickLabel',xtickLabels,'FontSize',8);
           %hT=rotateticklabel(gca,-90);
           h_legend=legend(handle,frsObj.fitResCell{frsObj.maxNumIndex}.lambda.dataLabels(fitNum),'Location','SouthEast');
           set(h_legend,'FontSize',10)
           hx=get(gca,'XLabel');  
           set([hx hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
           v=axis;
           frsObj.setCoeffRange(v(3),v(4));
           grid on;
           axis tight;
           if(frsObj.numCoeffs>1)
                xticklabel_rotate([],90,[],'Fontsize',10);%rotateticklabel(gca,-90); 
           end
            
        end
        function handle = plot3dCoeffSummary(frsObj,h)
            % handle = plot3dCoeffSummary(frsObj,h)
            % x-axis (covariate labels)
            % y-axis (coefficient values)
            % z-axis (histogram of coefficient values)
            if(nargin<2)
                h=gca;
            end
           frsObj.plotAllCoeffs(h);
           [N,edges] = frsObj.binCoeffs;
           numPlotCoeffs = length(frsObj.plotParams.xLabels);
           handle=ribbon(repmat(edges,[1 numPlotCoeffs]),N); 
           set(handle,'edgecolor','none');
           alpha(.6);
           legend off;
           view(gca,[71.5 28]);
           set(gca,'Box','off', 'Projection','perspective','Color',[0.831372549019608 0.815686274509804 0.784313725490196]);
           grid on; axis tight;
        end
        function handle = plot2dCoeffSummary(frsObj,h)
             % handle = plot2dCoeffSummary(frsObj,h)
             % histogram of regression coefficients for each unique
             % covariate.
            if(nargin<2)
                h=gca;
            end
            [N,edges,percentSig] = frsObj.binCoeffs;
            offset=0;
            numPlotCoeffs = length(frsObj.plotParams.xLabels);
            for i=1:numPlotCoeffs
                offset=offset+1;
                handle(i)=plot(h,edges,N(:,i)+offset); hold on;
%                 plot(edges,N(:,i)+offset,'.');
                 set(gca,'xtick',[],'ytick',[]);
            end
            % Reset the bottom subplot to have xticks
            set(gca,'xtickMode', 'auto');
            set(gca,'ytick',1:length(frsObj.plotParams.xLabels),'ytickLabel',frsObj.plotParams.xLabels,'FontSize',6);

            offset=0;
            for i=1:numPlotCoeffs
                offset=offset+1;
                text(frsObj.coeffRange.maxVal,offset,strcat(num2str(percentSig(i)*100,'%2.f'),'%_{sig}'),'FontSize',6); hold on;
            end
            
        end
        function handle = plotKSSummary(frsObj,neurons)
            % handle = plotKSSummary(frsObj,neurons)
            % For all of the distinct neurons in the the FitResSummary,
            % plots the corresponding KS plot with all the differentt fits
            % available for each neuron. Each neuron may have different numbers
            % of model fits: eg. Neuron j may have 3 fits while
            % neuron i may have 6.
            if(nargin<2||isempty(neurons))
               neurons = 1:frsObj.numNeurons;
            end
            if(max(neurons)>frsObj.numNeurons)
               error('Indices must be <= numNeurons'); 
            end
                
            handle = figure;
            numToPlot=length(neurons);
            cnt=0;
            for i=neurons
                cnt=cnt+1;
                if(numToPlot==1)
                    %dont subplot
                elseif(numToPlot<=2)
                    subplot(1,2,cnt);
                elseif(numToPlot<=4)
                    subplot(2,2,cnt)
                elseif(numToPlot<=8)
                    subplot(2,4,cnt)
                elseif(numToPlot<=12)
                    subplot(3,4,cnt)
                elseif(numToPlot<=16)
                    subplot(4,4,cnt)
                elseif(numToPlot<=20)
                    subplot(5,4,cnt)
                elseif(numToPlot<=24)
                    subplot(6,4,cnt)
                elseif(numToPlot<=40)
                    subplot(10,4,cnt)
                else
                    subplot(10,10,cnt)
                end
                frsObj.fitResCell{i}.KSPlot; 
                if(i~=neurons(end))
                    legend off; ylabel(''); xlabel(''); title('');
                else
                     h= legend; set(h,'Location','Best'); ylabel(''); xlabel(''); title('');
                end
                text(.4,.9,['N' num2str(i)]);
                set(gca,'xtick',[0 1],'ytick',[0 1])
            end
            
        end
        function handle = plotAIC(frsObj)
            % handle = plotAIC(frsObj)
            % Plot mean +/- 1 standard error from the mean for the AIC for
            % each fit.
           AICdata=frsObj.AIC;
           mData=nanmean(AICdata,1);
           numNeurons = size(AICdata,1);
           sData=nanstd(AICdata,0,1)./sqrt(numNeurons);
           ciUp = mData+sData;
           ciDown = mData-sData;
           
           x=1:frsObj.numResults;
           plot(x,mData,'r','Linewidth',3); hold all;
           faceColor='r';
           p=patch([x, fliplr(x)],[ciUp fliplr(ciDown)],faceColor);
           set(p,'facecolor',faceColor,'edgecolor','none');
           alpha(.5);     
%            set(gca,'xticklabelmode','auto','xtickmode','auto');
           set(gca,'xtick',x,'xticklabel',frsObj.fitNames);
           if(length(x)>1)
                xticklabel_rotate([],90,[],'Fontsize',8);%rotateticklabel(gca,-90); 
           end
           
        end
        
        function handle = plotBIC(frsObj)
            % handle = plotBIC(frsObj)
            % Plot mean +/- 1 standard error from the mean for the BIC for
            % each fit.
           BICdata=frsObj.BIC;
           mData=nanmean(BICdata,1);
           numNeurons = size(BICdata,1);
           sData=nanstd(BICdata,0,1)./sqrt(numNeurons);
           ciUp = mData+sData;
           ciDown = mData-sData;
           
           x=1:frsObj.numResults;
           plot(x,mData,'r','Linewidth',3); hold all;
           faceColor='r';
           p=patch([x, fliplr(x)],[ciUp fliplr(ciDown)],faceColor);
           set(p,'facecolor',faceColor,'edgecolor','none');
           alpha(.5);     
%            set(gca,'xticklabelmode','auto','xtickmode','auto');
           set(gca,'xtick',x,'xticklabel',frsObj.fitNames);
           if(length(x)>1)
                xticklabel_rotate([],90,[],'Fontsize',8);%rotateticklabel(gca,-90); 
           end
        end
        
        
        function handle = plotResidualSummary(frsObj)
            handle = figure;
            for i=1:frsObj.numNeurons
                if(frsObj.numNeurons<=4)
                    subplot(2,2,i)
                elseif(frsObj.numNeurons<=8)
                    subplot(2,4,i)
                elseif(frsObj.numNeurons<=12)
                    subplot(3,4,i)
                elseif(frsObj.numNeurons<=16)
                    subplot(4,4,i)
                elseif(frsObj.numNeurons<=20)
                    subplot(5,4,i)
                elseif(frsObj.numNeurons<=24)
                    subplot(6,4,i)
                elseif(frsObj.numNeurons<=40)
                    subplot(10,4,i)
                else
                    subplot(10,10,i)
                end
                frsObj.fitResCell{i}.plotResidual; 
                if(i~=frsObj.numNeurons)
                    legend off; ylabel(''); xlabel(''); title('');
                else
                    h= legend; set(h,'Location','BestOutside'); ylabel(''); xlabel(''); title('');
                end
            end
        end
            
        function handle = plotSummary(frsObj)
            % handle = plotSummary(frsObj)
            % 
            scrsz = get(0,'ScreenSize');
            handle=figure('OuterPosition',[scrsz(3)*.1 scrsz(4)*.1 scrsz(3)*.9 scrsz(4)*.9]);
            h1=subplot(2,4,[1 2 5 6]);frsObj.plotAllCoeffs(h1); grid off;
            title({'GLM Coefficients Across Neurons';'with 95% CIs (* p<0.05)'},'FontWeight','bold','FontSize',11,'FontName','Arial');
            
            %subplot(2,4,[2 3]);frsObj.plot3dCoeffSummary; %rotateticklabel(get(gca,'ytickLabels'),0);
            subplot(2,4,[3 4]); boxplot(frsObj.KSStats,frsObj.fitNames,'labelorientation','inline'); 
            ylabel('KS Statistics');               
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx hy],'FontName', 'Arial','FontSize',11,'FontWeight','bold');
            title('KS Statistics Across Neurons','FontWeight','bold','FontSize',11,'FontName','Arial');
%             subplot(2,4,[6 7]);frsObj.plot2dCoeffSummary; %rotateticklabel(get(gca,'ytickLabels'),0);
            
            subplot(2,4,7); frsObj.getDiffAIC(1); 
            ylabel('\Delta AIC'); %xticklabel_rotate([],45,[],'Fontsize',6);
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx hy],'FontName', 'Arial','FontSize',11,'FontWeight','bold');
            title('Change in AIC Across Neurons','FontWeight','bold','FontSize',11,'FontName','Arial');
            
            subplot(2,4,8); frsObj.getDiffBIC(1); 
            ylabel('\Delta BIC'); %xticklabel_rotate([],45,[],'Fontsize',6);
            hx=get(gca,'XLabel');  hy=get(gca,'YLabel');
            set([hx hy],'FontName', 'Arial','FontSize',11,'FontWeight','bold');
           

            title('Change in BIC Across Neurons','FontWeight','bold','FontSize',11,'FontName','Arial');
            
        end     
        function handle = boxPlot(frsObj,X,diffIndex,h,dataLabels,varargin)
            if(nargin<3)
                h=gca;
            end
            if(nargin<5 || isempty(dataLabels))
                [~,columns] = size(X);
                tempIndex = 1:frsObj.numResults;
                actIndex  = find(tempIndex~=diffIndex);
                
                if(~isempty(actIndex))
                    for i=1:columns
                        if(length(actIndex)==columns)
                            dataLabels{i} = [frsObj.fitNames{actIndex(i)} ' - ' frsObj.fitNames{diffIndex}];
                        end
                    end
                else
                    dataLabels{1}=frsObj.fitNames{diffIndex}; %only put the name of the fit since no other fits
                end
                   
                
                if(columns>1)
                    handle = boxplot(X,strvcat(dataLabels));
                else
                    handle = boxplot(X,dataLabels); % when only one column
                end
                
                
               % set(gca,'xticklabelmode','auto','xtickmode','auto');
               % set(gca,'xtick',1:length(dataLabels),'xticklabel',dataLabels);
               % hT=rotateticklabel(gca,90);
              %  FitResSummary.xticklabel_rotate([],45,[],'interpreter','none');
            elseif(nargin>5)
                handle = boxplot(h,X,strvcat(dataLabels),varargin{:});
            elseif(nargin==5)
                handle = boxplot(h,X,strvcat(dataLabels));
            end
            h = get(get(gca,'child'),'child');
            group_name_handle = findobj(h,'type','text');
            group_name_handle = flipud(group_name_handle); %place in correct order - find obj returns backwards
            v=axis;
            vdiffy = v(4)-v(3);
            vdiffx = v(2)-v(1);
            for j=1:length(group_name_handle)
                text(0,0,get(group_name_handle(j),'string'),'color','k','position',[j-.0*vdiffx v(3)-.02*vdiffy 0],'rotation',-90,'Fontsize',8);
            end
            delete(group_name_handle);
            
        end
        
        function structure = toStructure(frsObj)
            fNames = fieldnames(frsObj);
            for i=1:length(fNames)
               currObj = frsObj.(fNames{i});
               if(isa(currObj,'double')||isa(currObj,'cell'))
                   if(strcmp(fNames{i},'fitResCell'))
                       structure.(fNames{i}) = FitResult.CellArrayToStructure(frsObj.(fNames{i}));
                   else
                       structure.(fNames{i}) = frsObj.(fNames{i});
                   end
               end
                
            end
            
        end
        
        function [coeffIndex, epochId,numEpochs] = getCoeffIndex(frsObj,fitNum,sortByEpoch)
          if(nargin<3 || isempty(sortByEpoch))
             sortByEpoch=0; 
          end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:frsObj.numResults;
          end
          if(isempty(frsObj.plotParams))
               frsObj.computePlotParams;
          end   
          [histIndex, epochId] = frsObj.getHistIndex(fitNum,sortByEpoch);
          allIndex = 1:length(frsObj.uniqueCovLabels);
          
          nonHistIndex = setdiff(allIndex,histIndex);
          nonNANIndex = find(sum(~isnan(squeeze(frsObj.plotParams.bAct(:,fitNum,:))),2)>=1);
          actCoeffIndex = nonHistIndex(ismember(nonHistIndex, nonNANIndex));
          allCoeffTerms = frsObj.uniqueCovLabels(actCoeffIndex);
          epochStartInd=regexp(allCoeffTerms,'_\{\d*\}','start'); 
          epochEndInd=regexp(allCoeffTerms,'_\{\d*\}','end');
         
          allCoeffIndex = [];
          epochsExist =0;
          nonEpochIndex=[];
          for i=1:length(allCoeffTerms)
              if(~isempty(allCoeffTerms{i}))
                allCoeffIndex = [allCoeffIndex i];
                 if(~isempty(epochStartInd{i}))
                     epochsExist=1;
                     actStart = epochStartInd{i}+2;
                     actEnd   = epochEndInd{i}-1;
                     numEpoch(i) = str2num(allCoeffTerms{i}(actStart:actEnd));
                 else
                 
                    nonEpochIndex = [nonEpochIndex i];
                    numEpoch(i) = 0; % make terms that only appear once part of epoch 0.
                    
                 end
              end
              
          end
          
          
          if(epochsExist && ~sortByEpoch)
              totalEpochs = unique(numEpoch);
              coeffIndex = [];
              if(nargout>1)
                  epochId=[];
              end
              for i=1:length(totalEpochs)
                 coeffIndex = [coeffIndex, find(numEpoch==totalEpochs(i))]; 
                 if(nargout>1)
                    epochId = [epochId, totalEpochs(i)*ones(size(find(numEpoch==totalEpochs(i))))];
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
          

%           nonNANIndex = find(sum(~isnan(frsObj.plotParams.bAct(:,fitNum)),2)>=1);
%           coeffIndex = coeffIndex(ismember(coeffIndex, nonNANIndex));
%             
          if(nargout>2)
              numEpochs = length(unique(epochId));
          end
          
        end
        
        function h=plotCoeffsWithoutHistory(frsObj,fitNum,sortByEpoch,plotSignificance)
           if(nargin<4 || isempty(plotSignificance))
               plotSignificance=1;
           end
           if(nargin<3 || isempty(sortByEpoch))
              sortByEpoch = 0; 
           end
           if(nargin<2 || isempty(fitNum))
               fitNum = 1:frsObj.numResults;
           end
           if(isempty(frsObj.plotParams))
               frsObj.computePlotParams;
           end    
               
            
          coeffIndex = frsObj.getCoeffIndex(fitNum,sortByEpoch);
          h=frsObj.plotAllCoeffs([],fitNum,[],plotSignificance,coeffIndex);
          
          
            
        end
        
        function [histIndex, epochId,numEpochs] = getHistIndex(frsObj,fitNum,sortByEpoch)
            %if sortByEpoch==1 then we group all regression terms with the
            %same name one next to each other by epoch (time interval).
            %Otherwise, we show all epoch one terms, followed by all epoch
            %2 terms, etc.
           if(nargin<3 || isempty(sortByEpoch))
            sortByEpoch = 0;
           end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:frsObj.numResults;
          end
          if(isempty(frsObj.plotParams))
               frsObj.computePlotParams;
          end  
          
          
          allHistTerms = regexp(frsObj.uniqueCovLabels,'^[\w*');
          epochStartInd=regexp(frsObj.uniqueCovLabels,'\]_\{\d*\}','start'); 
          epochEndInd=regexp(frsObj.uniqueCovLabels,'\]_\{\d*\}','end');
          allHistIndex = [];
          epochsExist =0;
          for i=1:length(allHistTerms)
              if(~isempty(allHistTerms{i}))
                allHistIndex = [allHistIndex i];
                 if(~isempty(epochStartInd{i}))
                     epochsExist=1;
                     actStart = epochStartInd{i}+3;
                     actEnd   = epochEndInd{i}-1;
                     numEpoch(i) = str2num(frsObj.uniqueCovLabels{i}(actStart:actEnd));
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
             
                   
          nonNANIndex = find(sum(~isnan(frsObj.plotParams.bAct(:,fitNum)),2)>=1);
          histIndex = histIndex(ismember(histIndex, nonNANIndex));

          if(nargout>2)
              numEpochs = length(unique(epochId));
          end
            
        end
        function [coeffMat, labels, seMat] = getCoeffs(frsObj, fitNum)
             if(nargin<2 || isempty(fitNum))
                fitNum =1:frsObj.numResults;
            end
            sortByEpoch = 0; % Make sure we have different time series if the history is divided into epochs;
            [coeffIndex, epochId, numEpochs] = frsObj.getCoeffIndex(fitNum,sortByEpoch);
            epochNums = unique(epochId);
            
            
            coeffStrings = frsObj.uniqueCovLabels(coeffIndex);
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
               labels(i,epochIndices{i}) = frsObj.uniqueCovLabels(coeffStrIndex{i});
            end
            
            
%             for i=1:numEpochs
%                coeffIndMat(1:epochLength(i),i) = coeffIndex(epochIndices{i});
%                labels(1:epochLength(i),i) = frsObj.uniqueCovLabels(coeffIndMat(1:epochLength(i),i));
%             end
            
            
            
            if(length(fitNum)>1)
                coeffMat=cell(1,length(fitNum));
                seMat=cell(1,length(fitNum));
                for i=1:length(fitNum)
                    coeffMat{i} = nan(size(coeffIndMat,1),size(coeffIndMat,2), frsObj.numNeurons);
                    seMat{i} = nan(size(coeffIndMat,1),size(coeffIndMat,2), frsObj.numNeurons);
                    for j=1:length(uniqueCoeffs)
                        bTemp=squeeze(frsObj.plotParams.bAct(coeffStrIndex{j},i,:));  
                        seTemp=squeeze(frsObj.plotParams.seAct(coeffStrIndex{j},i,:));  
                        for k=1:frsObj.numNeurons
                            if(size(epochIndices,2)==1)
                                if(size(bTemp,2)==1)
                                    coeffMat(j,epochIndices{1},k) = bTemp(k);
                                    seMat(j,epochIndices{1},k) = seTemp(k);
                                else
                                    coeffMat(j,epochIndices{1},k) = bTemp(:,k);
                                    seMat(j,epochIndices{1},k) = seTemp(:,k);
                                end
                            else
                                if(size(bTemp,2)==1)
                                    coeffMat{i}(j,epochIndices{j},k) = bTemp(k);
                                    seMat{i}(j,epochIndices{j},k) = seTemp(k);
                                else
                                    coeffMat{i}(j,epochIndices{j},k) = bTemp(:,k);
                                    seMat{i}(j,epochIndices{j},k) = seTemp(:,k);
                                end
                            end
                        end
                    end
                end
            else
                coeffMat = nan(size(coeffIndMat,1),size(coeffIndMat,2),frsObj.numNeurons);
                seMat = nan(size(coeffIndMat,1),size(coeffIndMat,2),frsObj.numNeurons);
                for j=1:length(uniqueCoeffs)
                    bTemp=squeeze(frsObj.plotParams.bAct(coeffStrIndex{j},fitNum,:));   
                    seTemp=squeeze(frsObj.plotParams.seAct(coeffStrIndex{j},fitNum,:));   
                    for k=1:frsObj.numNeurons
                        if(size(epochIndices,2)==1)
                            if(size(bTemp,2)==1)
                                if(numel(bTemp)==numel(epochIndices{1}) && frsObj.numNeurons==1)
                                    coeffMat(j,epochIndices{1},k) = bTemp(:);
                                    seMat(j,epochIndices{1},k) = seTemp(:);
                                else
                                    coeffMat(j,epochIndices{1},k) = bTemp(k);
                                    seMat(j,epochIndices{1},k) = seTemp(k);
                                end
                            else
                                coeffMat(j,epochIndices{1},k) = bTemp(:,k);
                                seMat(j,epochIndices{1},k) = seTemp(:,k);
                            end
                        else
                            if(size(bTemp,2)==1)
                                if(numel(bTemp)==numel(epochIndices{j}) && frsObj.numNeurons==1)
                                    coeffMat(j,epochIndices{j},k) = bTemp(:);
                                    seMat(j,epochIndices{j},k) = seTemp(:);
                                else
                                    coeffMat(j,epochIndices{j},k) = bTemp(k);
                                    seMat(j,epochIndices{j},k) = seTemp(k);
                                end
                                
                            else
                                coeffMat(j,epochIndices{j},k) = bTemp(:,k);
                                seMat(j,epochIndices{j},k) = seTemp(:,k);
                            end
                        end
                    end
                end

            end
            if(frsObj.numNeurons==1)
                coeffMat=coeffMat';
                seMat=seMat';
            end
        end
        
        function [histMat, labels] = getHistCoeffs(frsObj,fitNum)
            if(nargin<2 || isempty(fitNum))
                fitNum =1:frsObj.numResults;
            end
            sortByEpoch = 0; % Make sure we have different time series if the history is divided into epochs;
            [histIndex, epochId, numEpochs] = frsObj.getHistIndex(fitNum,sortByEpoch);
            epochNums = unique(epochId);
            
            
            histcoeffStrings = frsObj.uniqueCovLabels(histIndex);
            baseStringEndIndex =regexp(histcoeffStrings,'_\{\d\}','start');
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
            for i=1:length(uniqueCoeffs)
               histcoeffIndMat(i,epochIndices{i}) = histcoeffStrIndex{i};
               labels(i,epochIndices{i}) = frsObj.uniqueCovLabels(histcoeffStrIndex{i});
            end
            
            
%             for i=1:numEpochs
%                coeffIndMat(1:epochLength(i),i) = coeffIndex(epochIndices{i});
%                labels(1:epochLength(i),i) = frsObj.uniqueCovLabels(coeffIndMat(1:epochLength(i),i));
%             end
            
            
            
            if(length(fitNum)>1)
%                 histMat = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2), length(fitNum));
                histMat = cell(1, length(fitNum));
                for i=fitNum
                    histMat{i} = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2), frsObj.numNeurons);
                    for j=1:length(uniqueCoeffs)
                        bTemp=squeeze(frsObj.plotParams.bAct(histcoeffStrIndex{j},i,:));   
                        for k=1:frsObj.numNeurons
                            histMat{i}(j,epochIndices{j},k) = bTemp(k,:);
                        end
                    end
                end
            else
                histMat = nan(size(histcoeffIndMat,1),size(histcoeffIndMat,2),frsObj.numNeurons);
                
                for j=1:length(uniqueCoeffs)
                    bTemp=squeeze(frsObj.plotParams.bAct(histcoeffStrIndex{j},fitNum,:));
                    for k=1:frsObj.numNeurons
                        histMat(j,epochIndices{j},k) = bTemp(k,:);
                    end
                end

            end
            
            

            
        end
        function h=plotHistCoeffs(frsObj,fitNum,sortByEpoch,plotSignificance)
          if(nargin<4 || isempty(plotSignificance))
            plotSignificance=1;
          end
          if(nargin<3 || isempty(sortByEpoch))
             sortByEpoch=0; 
          end
          if(nargin<2 || isempty(fitNum))
              fitNum = 1:frsObj.numResults;
          end
          if(isempty(frsObj.plotParams))
               frsObj.computePlotParams;
          end  
          histIndex = frsObj.getHistIndex(fitNum,sortByEpoch);
          h=frsObj.plotAllCoeffs([],fitNum,[],plotSignificance,histIndex);
        end
        
        
        
        
    end 
    
    methods (Static)
        function frsObj = fromStructure(structure)
            resultsCell = FitResult.fromStructure(structure.fitResCell);
            frsObj = FitResSummary(resultsCell);
        end
        
    end
    
end

%% Helper Functions
function dMat = computeDiffMat(MAT,diffIndex)

     columns = size(MAT,2);
     index = 1:columns;
     dMat = MAT(:,index~=diffIndex) - MAT(:,diffIndex)*(ones(1,columns-length(diffIndex)));
            
end
function [uniqueLabels, indexIntoOriginal, restoreIndex] = getUniqueLabels(covLabels)
            offset = 0;
            for i=1:length(covLabels)
                currLabels = covLabels{i};                
                allLabels((1:length(currLabels))+offset) = currLabels;
                offset=length(allLabels);
            end
            [uniqueLabels, indexIntoOriginal, restoreIndex] = unique(allLabels);
end
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
    if (nargin < 3 || isempty(varargin{1})) & (~exist('XTick') | isempty(XTick)),
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
    textsizes = cell2mat(get(hText,'extent'))       ;
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

