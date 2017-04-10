classdef TrialConfig <handle
%TRIALCONFIG holds the configuration for a particular analysis
% Usage:
% >>tcObj=TrialConfig(covMask,sampleRate, history,ensCovHist,covLag,name)
%
% All parameters are optional
% <a href="matlab: methods('TrialConfig')">methods</a>
% <a href="matlab:web('TrialConfigExamples.html', '-helpbrowser')">TrialConfig Examples</a> 
%
% see also <a href="matlab:help('Trial')">Trial</a>, 
% <a href="matlab:help('CovColl')">CovColl</a>, 
% <a href="matlab:help('nstColl')">nstColl</a>
%
% Reference page in Help browser
% <a href="matlab: doc('TrialConfig')">doc TrialConfig</a>

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
        covMask % which covariates are going to be used for analysis
        covLag  % deltaT for the covariate. Negative to shift back in time
        sampleRate % sampleRate to be used for the fitting of this trial
        history % History object this trial
        ensCovHist % History term for neiboring neurons - aka Ensemble History
        ensCovMask % Mask Matrix for neighboring neurons;
        name % Name of the configuration
        
    end
    
    
    methods
        function tcObj=TrialConfig(covMask,sampleRate, history,ensCovHist,ensCovMask,covLag,name)
           % tcObj=TrialConfig(covMask,sampleRate, history,covHist,covLag,name)
           % covMask: List of covariate names and labels that are to be used:
           %  ex: {{'Position','x'},{'Velocity','v_x'}} 
           %        would use the component x of the Position covariate and
           %        the component v_x of the Velocity covariate for the fits
           %        - if empty then no covariates are used
           %
           % sampleRate: if empty uses the current sampleRate. Otherwise
           %             resamples all of the data to the specified sample
           %             rate.
           %
           % history: history object, cell of history objects,  or a vector of windowTimes
           % ensCovHist: same as history. Used to determine ensemble covariate history
           % covLag : delay to be used with covariates
           % name:    name o f the configuration
           if(nargin<7 || isempty(name))
               name ='';
           end
           if(nargin<6 || isempty(covLag))
                covLag=[];
           end
           if(nargin<5 || isempty(ensCovMask))
               ensCovMask=[];
           end
           if(nargin<4 || isempty(ensCovHist))
               ensCovHist = [];
           end
           if(nargin<3 || isempty(history))
               history = [];
           end
           if(nargin<2 || isempty(sampleRate))
               sampleRate =[];
           end
           if(nargin<1 || isempty(covMask))
               covMask=[];
           end
           
           tcObj.covMask = covMask;
           tcObj.sampleRate = sampleRate;
           tcObj.history = history;
           tcObj.ensCovHist = ensCovHist;
           tcObj.ensCovMask = ensCovMask;
           tcObj.covLag  = covLag;
           tcObj.name    = name;
        end
        function setConfig(tcObj, trial)
        % setConfig(tcObj, trial)
        % applies this configuration to trial.
            if(~isempty(tcObj.history))
                trial.setHistory(tcObj.history);
            else
                trial.resetHistory;
            end

        
            if(~isempty(tcObj.sampleRate))
                if(trial.sampleRate~=tcObj.sampleRate)
                    trial.resample(tcObj.sampleRate);
                end
            else
                %trial.restoreToOriginal;
                %keep sampleRate
            end
            
            %if(~isempty(tcObj.covMask))
            trial.setCovMask(tcObj.covMask);
                %all mask to be empty
            %else
            %    trial.resetCovMask;
            %end
            
            if(~isempty(tcObj.covLag))
                trial.shiftCovariates(tcObj.covLag);
            end
            
            
            if(~isempty(tcObj.ensCovHist))
                trial.setEnsCovHist(tcObj.ensCovHist);
                trial.setEnsCovMask(tcObj.ensCovMask);
            else
                trial.setEnsCovHist; %sets it to be empty
                trial.setEnsCovMask; %sets it to the default;
            end
%             trial.setTrialTimesFor('training');
        end
        function n=getName(tcObj)
            % n=getName(tcObj)
            % returns a string with the name of this configuration
            n = tcObj.name;
        end
        function setName(tcObj,n)
            % setName(tcObj,n)
            % sets the name of the current configuration to n
            tcObj.name = n;
        end
        function structure = toStructure(tcObj)
            % structure = toStructure(tcObj)
            % structure is a standard matlab structure that can be saved as
            % a .mat file. 
            % A TrialConfig object can be obtained by calling 
            % tcObj = TrialConfig.fromStructure(structure);
            fNames = fieldnames(tcObj);
            for i=1:length(fNames)
               structure.(fNames{i}) = tcObj.(fNames{i}); 
            end
            
        end
        
    end
    methods (Static)
        function tcObj = fromStructure(structure)
            % tcObj is a TrialConfig object that is reconstructed from the
            % structure
            tcObj=TrialConfig(structure.covMask,structure.sampleRate, ...
                              structure.history,structure.ensCovHist, ...
                              structure.covLag,structure.name);
        end
    end
    
    
    
end
