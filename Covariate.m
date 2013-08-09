classdef Covariate < SignalObj
% COVARIATE Covariates are signals (of class SignalObj) with a mean mu and a standard
% deviation sigma.
%
% cov  = Covariate(time, data, name, xlabelval, xunits, yunits, dataLabels, plotProps)
% All inputs are passed to the superclass SignalObj.
%
% Each dimenion of a covariate signal has a mean. 
% cov.mu    - SignalObj reprenting mean of each component over time
% cov.sigma - SignalObj reprenting standard deviation of each component
%             over time.
%
% cov.getSigRep('standard') or cov.getSigRep is the original data 
% can also just use cov for the standard representation
% cov.getSigRep('zero-mean') is a zero mean version of the Signal 
%
% <a href="matlab: methods('Covariate')">methods</a>
% <a href="matlab:web('CovariateExamples.html', '-helpbrowser')">Covariate Examples</a> 
%
% see also <a href="matlab:help('SignalObj')">SignalObj</a>, <a href="matlab:help('CovColl')">CovColl</a>  
%
% Reference page in Help browser
% <a href="matlab: doc('Covariate')">doc Covariate</a>

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
    
    properties (Dependent = true)
          mu    %SignalObj representing the mean of each component of the Covariate across time
          sigma %Standard deviation of the covariate across time
          
    end
    properties
       ci    %a Confidence Interval object for the covariate 
    end
    
    methods
        function cov  = Covariate(varargin)
            %cov  = Covariate(time, data, name, xlabelval, xunits, yunits,
            %dataLabels, plotProps)
            cov@SignalObj(varargin{:});

        end  
        function newCov = computeMeanPlusCI(covObj,alphaVal)
            if(nargin<2)
                alphaVal=.05;
            end
             for k=1:length(covObj.time)
              [f,x] = ecdf(squeeze(covObj.data(k,:)));
              CIs(k,1) = x(find(f<alphaVal/2,1,'last'));
              CIs(k,2) = x(find(f>(1-alphaVal/2),1,'first'));
             end
            confInt = ConfidenceInterval(covObj.time,CIs,'CI','time','s','');
            newCov = mean(covObj,2);
            newCov.setConfInterval(confInt);
        end
        function h = plot(covObj,varargin)
           h=plot@SignalObj(covObj,varargin{:}); 
           if(covObj.isConfIntervalSet)
               handles = get(gca,'Children');
%                actHandles = [];
%                for k=1:length(handles);
%                    if(strcmp(get(handles(k),'type'),'line'))%&& ~isempty(get(handles(k),'DisplayName')))
%                        actHandles = [actHandles;handles(k)];
%                    end
%                end
               actHandles = handles(strcmp('line',get(get(gca,'Children'),'type'))); 
               s=get(actHandles);
              
               selectorArray = find(covObj.dataMask==1);
               for i=1:length(selectorArray)
                   actIndex = length(selectorArray)-(i-1);
                   TempColor=s(actIndex).Color;
                   covObj.ci{selectorArray(i)}.plot(TempColor);
                   
               end
               axis tight;
           end
            
        end
        function cov = getSubSignal(covObj,varargin)
            cov = getSubSignal@SignalObj(covObj,varargin);
            if(covObj.isConfIntervalSet)
                origIndex = zeros(1,cov.dimension);
                for i=1:cov.dimension
                    origIndex(i) = find(strcmp(cov.dataLabels{i},covObj.dataLabels));
                end
                cov.setConfInterval(covObj.ci(origIndex));
            end
        end
        function cSig = getSigRep(covObj,repType)
            % cSig = getSigRep(covObj,repType)
            % repType: 'standard'  - original representation
            %          'zero-mean' - zero mean representation
            if(nargin<2)
                repType = 'standard';
            end
            if(strcmp(repType, 'zero-mean'))
                cSig = covObj-covObj.mu;
            elseif(strcmp(repType,'standard'))
                cSig = covObj;
            else
                error('repType must be either ''zero-mean'' or ''standard'' ');
            end
        end
        function mu = get.mu(covObj)
            mu = mean(covObj);
        end
        function sigma = get.sigma(covObj)
            sigma = std(covObj);
        end
        
        function cov = filtfilt(covObj,varargin)
           cov=filtfilt@SignalObj(covObj,varargin{:});
        end
        
        
        function structure = toStructure(covObj)
            fNames = fieldnames(covObj);
            for i=1:length(fNames)
                currObj = covObj.(fNames{i});
                
                if(strcmp(fNames{i},'ci'))
                    if(covObj.isConfIntervalSet)
                        if(isa(covObj.ci,'ConfidenceInterval'))
                            structure.ci = covObj.ci.dataToStructure;
                        elseif(isa(covObj.ci,'cell'))
                          for j=1:length(covObj.ci)
                            ciTemp = covObj.ci{j};
                            structure.ci{j} = ciTemp.dataToStructure;
                          end
                        end
                    end
                elseif(isa(currObj,'double')||isa(currObj,'cell')||isa(currObj,'char'))
                    structure.(fNames{i}) = currObj;
                elseif(isa(currObj,'Covariate'))
                    structure.(fNames{i}) = currObj.dataToStructure;
                end
            end
            
              
            
        end
        
        function ans = isConfIntervalSet(covObj)
            ans = ~isempty(covObj.ci);
        end
        function setConfInterval(covObj, ciObj)
            if(isa(ciObj,'cell'))
                covObj.ci = ciObj;
            elseif(isa(ciObj,'ConfidenceInterval'))
                covObj.ci = {ciObj};
            end
        end
        function covOut = copySignal(covObj)
            covOut=copySignal@SignalObj(covObj);
            if(covObj.isConfIntervalSet)
                covOut.setConfInterval(covObj.ci);
            end
        end
        
        function covOut = plus(covObj,covObj2)
            covOut=plus@SignalObj(covObj,covObj2);
            
            if(isa(covObj,'Covariate')&&isa(covObj2,'Covariate'))
                if(covObj.isConfIntervalSet && ~covObj2.isConfIntervalSet)
                    for i=1:covObj.dimension
                        tempCi{i} = covObj.ci{i} + covObj2.getSubSignal(i);
                    end
                    covOut.setConfInterval(tempCi);
                elseif(covObj.isConfIntervalSet && covObj2.isConfIntervalSet)
                    for i=1:covObj.dimension
                        tempCi{i} = covObj.ci{i} + covObj2.ci{i};
                    end
                    covOut.setConfInterval(tempCi);
                elseif(~covObj.isConfIntervalSet && covObj2.isConfIntervalSet)
                    for i=1:covObj2.dimension
                        tempCi{i} = covObj2.ci{i} + covObj.getSubSignal(i);
                    end
                    covOut.setConfInterval(tempCi);
                end
            end
            
        end
        
        function covOut = minus(covObj,covObj2)
            covOut=minus@SignalObj(covObj,covObj2);
            
            if(isa(covObj,'Covariate')&&isa(covObj2,'Covariate'))
                if(covObj.isConfIntervalSet && ~covObj2.isConfIntervalSet)
                    for i=1:covObj.dimension
                        tempCi{i} = covObj.ci{i} - covObj2.getSubSignal(i);
                    end
                    covOut.setConfInterval(tempCi);
                elseif(covObj.isConfIntervalSet && covObj2.isConfIntervalSet)
                    for i=1:covObj.dimension
                        tempCi{i} = covObj.ci{i} - covObj2.ci{i};
                    end
                    covOut.setConfInterval(tempCi);
                elseif(~covObj.isConfIntervalSet && covObj2.isConfIntervalSet)
                    for i=1:covObj2.dimension
                        tempCi{i} = -covObj2.ci{i} + covObj.getSubSignal(i);
                    end
                    covOut.setConfInterval(tempCi);
                end
            end 
            
        end
        
        
        function covOut = dataToStructure(covObj)
            covOut=dataToStructure@SignalObj(covObj);
        end
        
    end
   methods (Static)
       function cov = fromStructure(structure)

            cov=Covariate(structure.time, structure.data, structure.name, structure.xlabelval, structure.xunits, structure.yunits, structure.dataLabels, structure.plotProps);
            fnames = fields(structure);
            if(any(strcmp('ci',fnames)))
                if(~isempty(structure.ci))
                    if(isa(structure.ci,'cell'))
                        for i=1:length(structure.ci)
                            ciTemp{i} = ConfidenceInterval.fromStructure(structure.ci{i});
                        end
                        cov.setConfInterval(ciTemp);
                    elseif(isa(structure.ci,'struct'))
                        cov.setConfInterval(ConfidenceInterval.fromStructure(structure.ci));
                    end
                end
            end
        
       end
       
   end
   
end

