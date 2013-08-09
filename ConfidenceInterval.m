classdef ConfidenceInterval < SignalObj
%ConfidenceInterval-represents the confidence interval for a time
%series or <a href="matlab:help('Covariate')">Covariate</a>.
% <a href="matlab: methods('ConfidenceInterval')">methods</a>
% Reference page in Help browser
% <a href="matlab: doc('ConfidenceInterval')">doc ConfidenceInterval</a>

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
        color
        value % percent of the CI .eg. 1-alpha - default 95%
    end
    
    methods 
        function ciObj = ConfidenceInterval(varargin)
            ciObj@SignalObj(varargin{:});
            ciObj.color = 'b';  
            ciObj.value = .95;
        end
        function setColor(ciObj,color)
           ciObj.color = color; 
            
        end
        function setValue(ciObj,value)
            ciObj.value = value;
        end
        
        function plot(ciObj,color,alphaVal,drawPatches)
            if(nargin<4)
                drawPatches=0;
            end
            if(nargin<3)
                alphaVal=.2;
            end
            if(nargin<2)
                color = ciObj.color;
            end
            ciData = ciObj.dataToMatrix;
            ciHigh=ciData(:,2);         
            ciLow=ciData(:,1);  
            time = ciObj.time;

            hold on;
            if(drawPatches==1)
                if(isa(color,'char'))
                    p=patch([reshape(time,[1 length(time)]) reshape(fliplr(time'),[1 length(time)])],...
                            [reshape(ciLow',[1 length(time)])       reshape(fliplr(ciHigh'),[1,length(time)])],strcat('',color,''));
                    set(p,'facecolor',color,'edgecolor','none');
                elseif(isa(color,'double'))
                    p=patch([reshape(time,[1 length(time)]) reshape(fliplr(time'),[1 length(time)])],...
                            [reshape(ciLow',[1 length(time)])       reshape(fliplr(ciHigh'),[1,length(time)])],color);
                    set(p,'facecolor',color,'edgecolor','none');                
                end
                alpha(alphaVal);
            else 
                if(isa(color,'char'))
                    p=plot(time, ciData);
                elseif(isa(color,'double'))
                    p=plot(time, ciData);
                    set(p,'Color',color); 

                end
                alpha(alphaVal);
            end
        end
    end
    methods (Static)
        function ciObj = fromStructure(structure)
            ciObj=ConfidenceInterval(structure.time, structure.signals.values, structure.name, structure.xlabelval, structure.xunits, structure.yunits, structure.dataLabels, structure.plotProps);
        end 
    end
end

