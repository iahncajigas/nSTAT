classdef Events
% EVENTS Events represent times of importance during an experiment that
% need to be highlighted in figures or used to identify certain epochs
% in the data. 
% 
% e = Events(eventTimes, eventLabels)
% e = Events(eventTimes, eventLabels,eventColor)
% eventTimes: a vector of times at which certain events occur.
% eventLabels: a cell array of strings containing the names of each
%              of the events indicated by eventTimes.
% eventColor: strings indicating the color for this event. Same as
%             colors strings in matlab's standard plot routine. If not
%             specified, default is red.
%
% The length of eventTimes and eventLabels must match for an Events
% object to be successfully created.
% <a href="matlab: methods('Events')">methods</a>
% <a href="matlab:web('EventsExamples.html', '-helpbrowser')">Events Examples</a> 
%
% Reference page in Help browser
% <a href="matlab: doc('Events')">doc Events</a>


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
        eventTimes; % time  of each event
        eventLabels;% label of each event
        eventColor; % color for plotting event
    end
    
    methods 
        function e = Events(eventTimes, eventLabels,eventColor)
            % e = Events(eventTimes, eventLabels,eventColor)
            if(nargin<3)
                eventColor = 'r';
            end
            
            if(length(eventTimes)==length(eventLabels))
                e.eventTimes  = eventTimes;
                e.eventLabels = eventLabels;
                e.eventColor  = eventColor;
            else
                error('Number of eventTimes must equal number of eventLabels');
            end
        end
        
        function h = plot(EObj,handle,colorString)
            % h = plot(EObj,handle)
            % plots the event on the current on the figure/subplot
            % specified by handle.
            % If handle is not specified, then handle defaults to the
            % current figure window.
            if((nargin<3) || isempty(colorString))
                colorString=EObj.eventColor;
            end
            if((nargin<2) || isempty(handle))
                handle=gca;
            end
            
            
            
            for j=1:length(handle) %An event can be told to plot itself on multiple axes
                set(gcf,'CurrentAxes',handle(j));
                v=axis;
%                 for k=1:length(EObj.eventTimes) 
%                     times   =   ones(1,4)*EObj.eventTimes(k) + 4*[-.001 .001 .001 -.001]*(v(2)-v(1));
%                     y       =   [v(3), v(3),v(4),v(4)];
                    times = repmat(EObj.eventTimes,[2 1]);
                    y =  [v(3); v(4)]*ones(1,length(EObj.eventTimes));
                    hold on;
                    h=plot(handle(j),times,y,'r','LineWidth',4);
                    %tcolor(1,1,1:3) = [1 0 0];
%                     p=patch(times,y,colorString);%'FaceAlpha',.7,'EdgeAlpha',0)
%                     set(p,'facecolor',colorString,'edgecolor','none');
%                     alpha(.5);
%                 end
                    % Create textbox
                v=axis;
                for i=1:length(EObj.eventTimes)
                    if( ((EObj.eventTimes(i)-v(1))/(v(2)-v(1))>=0) && (EObj.eventTimes(i)<=v(2)))
                        text((EObj.eventTimes(i)-v(1))/(v(2)-v(1))-.02,1.03,EObj.eventLabels{i},'rotation',0,'FontSize',10,'Color',[0 0 0],'Units', 'normalized'); %write event labels
                        %'Interpreter','latex'
                    end
                end
            end
        end
        
        function structure = toStructure(EObj)
           fNames = fieldnames(EObj);
           for i=1:length(fNames)
               structure.(fNames{i}) = EObj.(fNames{i});
           end
           
           
        end
        
            
    end
    methods (Static)
        function EObj = fromStructure(structure)
            if(~isempty(structure))
                fNames = fieldnames(structure);
                reqNames = {'eventTimes','eventLabels','eventColor'};
                for i=1:length(reqNames)
                   if(~any(strcmp(reqNames{i},fNames)))
                       error('Missing field in structure. Cant creats Events object!');
                   end
                end
                EObj=Events(structure.eventTimes,structure.eventLabels,structure.eventColor);
            else
                EObj = [];
            end
        end
    end
        
        
end
    
function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from axis to figure coords
% Transforms [axx axy] or [xypos] from axes hAx (data) coords into coords
% wrt GCF for placing annotation objects that use figure coords into data
% space. The annotation objects this can be used for are
%    arrow, doublearrow, textarrow
%    ellipses (coordinates must be transformed to [x, y, width, height])
% Note that line, text, and rectangle anno objects already are placed
% on a plot using axes coordinates and must be located within an axes.
% Usage: Compute a position and apply to an annotation, e.g.,
%   [axx axy] = ginput(2);
%   [figx figy] = getaxannopos(gca, axx, axy);
%   har = annotation('textarrow',figx,figy);
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')'])

%% Obtain arguments (only limited argument checking is performed).
% Determine if axes handle is specified
if length(varargin{1})== 1 && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end);
else
	hAx = gca;
end;
% Parse either a position vector or two 2-D point tuples
if length(varargin)==1	% Must be a 4-element POS vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Two tuples (start & end points)
end
%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));
%% Transform data from figure space to data space
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
	varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
	pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
	pos(3) = pos(3)*axpos(3)/axwidth;
	pos(4) = pos(4)*axpos(4)/axheight;
	varargout{1} = pos;
end
%% Restore axes units
set(hAx,'Units',axun)
end



    