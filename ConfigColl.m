classdef ConfigColl < handle
%CONFIGCOLL - a collection of different configurations
% <a href="matlab: methods('ConfigColl')">methods</a>
% <a href="matlab:web('ConfigCollExamples.html', '-helpbrowser')">ConfigColl Examples</a> 
%
% see also <a href="matlab:help('TrialConfig')">TrialConfig</a>, <a
% href="matlab:help('CovColl')">Trial</a>, <a
% href="matlab:help('Analysis')">Analysis</a>
%
% Reference page in Help browser
% <a href="matlab: doc('ConfigColl')">doc ConfigColl</a>
    

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
        numConfigs % number of configurations in this collection
        configNames% names of the configurations 
    end

    properties (Hidden) 
        configArray
    end
    
    methods
       function tcColl = ConfigColl(tcObj)
            % tcColl = ConfigColl(tcObj)
            % creates a ConfigColl object from a TrialConfig objects
            % tcObj can be a cell of TrialConfigs, a cell of strings (just
            % names)
            if(nargin<1)
                tcObj = [];
            end
             tcColl.numConfigs = 0;
             tcColl.configArray=[];
             tcColl.addConfig(tcObj);
       end        
       function addConfig(tcColl,tcObj)
           if(~isempty(tcObj))
               if(isa(tcObj,'cell'))
                   if(isa(tcObj{1},'TrialConfig'))
                       for i=1:length(tcObj)
                           tcColl.addConfig(tcObj{i});
                       end
                   elseif(isa(tcObj{1},'char'))
                       for i=1:length(tcObj)
                           tcColl.addConfig(tcObj{i});
                       end
                   elseif(isempty(tcObj{1}))
                       for i=1:length(tcObj)
                           tcColl.addConfig(tcObj{i});
                       end
                   end
                   
               elseif(isa(tcObj,'TrialConfig'))
                    tcColl.numConfigs = tcColl.numConfigs+1;
                    tcColl.setConfigNames(tcObj.name,tcColl.numConfigs);
                    tcColl.configArray{tcColl.numConfigs} = tcObj;
               elseif(isa(tcObj,'char'))
                    tcColl.numConfigs = tcColl.numConfigs+1;
                    tcColl.setConfigNames(tcObj.name,tcColl.numConfigs);
                    tcColl.configArray{tcColl.numConfigs} = tcObj;
               end
           else
               tcColl.numConfigs = tcColl.numConfigs+1;
               tcColl.setConfigNames('Empty Config',tcColl.numConfigs);
               tcColl.configArray{tcColl.numConfigs} = {'Empty Config'};
           end
       end            
       function config = getConfig(tcColl,index)
            %config = getConfig(tcColl,index)
            % retunrs the TrialConfig specified by index only if index is
            % within 0 to numConfigs. Otherwise returns an error 'Index out
            % of bounds'
            if(index>0 && index<=tcColl.numConfigs)
                config = tcColl.configArray{index};
            else
                error('Index Out of Bounds');
            end
       end       
       function setConfig(tcColl,trial,index)
           % setConfig(tcColl,trial,index)
           % Set the configuration specified by the index to Trial
            config=tcColl.getConfig(index);
            if(isa(config,'TrialConfig'))
                config.setConfig(trial);
            else
                error('Cannot Set Empty Configs');
            end
       end
       function cArray = getConfigNames(tcColl,index)
           % cArray = getConfigNames(tcColl,index)
           % returns a cell array of strings with the names of the
           % TrialConfigs specified by index.
           if(nargin<2)
               index=1:tcColl.numConfigs;
           end
           cArray=cell(1,length(index));
           for i=1:length(cArray)
               tempName=tcColl.configNames{index(i)};
               if(isempty(tempName))
                   cArray{i} = ['Fit ' num2str(i)];
               else
                   cArray{i} = tempName;
               end
             
           end
           
       end
       function setConfigNames(tcColl, names, index)
           % setConfigNames(tcColl, names, index)
           % sets the TrialConfigs specified by index to have the names
           % specified in names.
           % if names is a string, then index must be length 1
           % if names is a cell array with n string entries, then index
           % must be length n.
           if(nargin<3)
               index=1:tcColl.numConfigs;
           end
           if((max(index)<=tcColl.numConfigs) && min(index)>0)
               if(isa(names,'char') && length(index)==1)
                   if(isempty(names))
                        tcColl.configNames{index} = ['Fit ' num2str(tcColl.numConfigs)];
                   else
                        tcColl.configNames{index} = names;
                   end
               elseif(isa(names,'cell') && isa(names{1},'char'))
                   if(length(index)==length(names))
                       for i=1:length(names)
                           tcColl.configNames{index(i)} = names{i};
                       end
                   else
                       error('If specifying multiple names, either names must be the same size as numConfigs or and index specified');
                   end
               end
           end
       end
       function subsetConfigs=getSubsetConfigs(tcColl, subset)
           for i=1:length(subset)
               tempconfigs{i} = tcColl.getConfig(subset(i));
           end
           subsetConfigs=ConfigColl(tempconfigs);
       end
       function structure = toStructure(tcColl)
          fNames = fieldnames(tcColl); 
          for i=1:length(fNames)
             structure.(fNames{i}) = tcColl.(fNames{i});
          end
          for i=1:tcColl.numConfigs
             structure.configArray{i} = tcColl.configArray{i}.toStructure; 
          end
          
       end

    end
    methods (Static)
        function tcColl = fromStructure(structure)
            c=cell(1,structure.numConfigs);
            for i=1:structure.numConfigs
               c{i} = TrialConfig.fromStructure(structure.configArray{i});
            end
            tcColl = ConfigColl(c);
        end
        
    end
    
end
