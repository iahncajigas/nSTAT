classdef SignalObj < handle
%SIGNALOBJ Class representing a signal abstraction
%   SignalObj consist of data that is indexed by time (as a default). The
%   indexing variable can be any other type of data and the x-axis labels 
%   modified to represent this change. 
%
%   A SignalObj can be multivariate in that the data can have more than one component. The
%   sample rate of the SignalObj is determined by the time increment used
%   in the time sequence used when the SingalObj is created
%
%   Usage:  
%   >> s=SignalObj(time, data, name, xlabelval, xunits, yunits, dataLabels, plotProps)
%
%   Only time and data need to be specified. Other arguments are optional.
%   
%
%  time: indexing variable for the data. n x 1 or 1 x n array. The sample 
%        rate is determined by the time increment between samples of this vector. Units of 
%        [sec] are assumed, but need not be used. If the time vector is 
%        in units of [sec], the sampleRate is in units of [Hz]. If the 
%        time vector is in units of [msec] then the sampleRate is in 
%        units of [1/msec] or 10^3 [Hz].
%
%  data: n x m or m x n array reprenting the signal at each index of the time vector.
%        The dimension that is compatible with the time vector will be automatically detected. 
%        Thus a SignalObj can be created by either passing the data matrix or its transpose. The remaining
%        dimension will determine the dimensionality of the SignalObj.
%
%  name: string that determines the name of the signal. This is used to
%        label the y-axis of the SignalObj.
%
%  xlabelval: A string specifying the name of the indexing variable. If
%        this value is not specified, 'time' is used.
%
%  xunits: A string specifying the name of the units of the indexing
%          variable. In not specified, 'sec' is used.
%         
%  yunits: A string specifying the units of the SignalObj. Used when plotting the SignalObj.       
%
%  dataLabels: If data is multivariate, the names of the components of the SignalObj can be specified. 
%              These can be used to reference specific data within the
%              signal (e.g. the x-component of a 3-d vector) and are
%              also used for plotting. SignalObj's will be created for
%              each component of the orignal SignalObj under the
%              vars field. Can be specified all at once or by a cell of
%              strings.
%
%  plotProps:  Can be specified for each component of the SignalObj
%              individually or by a cell of string of same dimension as the 
%              number of components in the data. 
%
%
% <a href="matlab: methods('SignalObj')">methods</a>
% <a href="matlab:web('SignalObjExamples.html', '-helpbrowser')">SignalObj Examples</a> 
%
% Reference page in Help browser
% <a href="matlab:doc('SignalObj')">doc SignalObj</a>
    

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
        name       % name of the SignalObj
        time       % time vector. Time increment determines sampleRate
        data       % actual SignalObj data
        dimension  % number of different components of the SignalObj
        minTime    % minimum Time value of the SignalObj
        maxTime    % maximum Time value of the SignalObj
        xlabelval  % label to use for the x-axis
        xunits     % units for x-axis
        yunits     % units for y-axis data
        dataLabels % labels for each dimension of the data;
        dataMask   % vector same length as SignalObj dimension. a 1 indicates this SignalObj should be output a 0 otherwie
        sampleRate % Hz if time is in seconds
        plotProps  % Plotting properties        
    end
    properties (Hidden)
        origSampleRate
        originalTime %original timeVector
        originalData %original Data
    end
%     properties (Dependent = true)    
%         vars %Contains subfields of the same names as the dataLabels that contain Signals with only the data corresponding to that label.
%     end
    methods
        %Constructor
         function s=SignalObj(time, data, name, xlabelval, xunits, yunits, dataLabels, plotProps)

            if(nargin<6)
                yunits='';
            end
            if(nargin<5)
                xunits='s';
            end
            if(nargin<4)
                xlabelval='time';
            end
            if(nargin<3)
                name='';
            end
            [l,w]=size(time);
            if(l>=w);
                if(w>1)
                    error('Time vector can only have one dimension');
                else
                    s.time=time;
                    
                end
            elseif(l<=w)
                if(l>1)
                    error('Time vector can only have one dimension');
                else
                    s.time=time';
                end
            end
            s.originalTime=s.time;
            [l,w]=size(data);
           if(l==length(s.time));
                s.data=data;
                s.dimension =w;
            elseif(w==length(s.time));
                s.data=data';
                s.dimension=l;
           else
               error('Data dimensions do not match the time vector specified');
           end
           s.originalData = s.data;
            if(nargin <7)
                if(s.dimension==0)
                    dataLabels ='';
                else
                    for i=1:s.dimension
                        dataLabels{i} = '';
                    end
                end
            end
            s.dataMask  = ones(1,s.dimension);   
            if(nargin<8)
                plotProps = cell(s.dimension,1);
            end
            deltaT = round(1000*mean(diff(s.time)))/1000; %To avoid round-off error, when computing samplerate
            s.sampleRate = 1/deltaT; 
            s.origSampleRate = s.sampleRate;
            s.name=name;
            s.xlabelval=xlabelval;
            s.xunits=xunits;
            s.yunits=yunits;
            s.minTime=min(s.time);
            s.maxTime=max(s.time);
            s.setPlotProps(plotProps);
            s.setDataLabels(dataLabels);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set functions    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setName(sObj,name) 
            % setName(sObj,name) 
            % set the name after construction
            if(isa(name,'char'))
                sObj.name = name;
            else
                error('Name must be a string!');
            end
        end     
        function setXlabel(sObj,name)
            %setXlabel(sObj,name)
            %set the x-axis label to string name
            sObj.xlabelval = name;
        end
        function setYLabel(sObj,name)
            %setYLabel(sObj,name)
            %set the ylabel to string name;
            %Same as calling setName(sObj,name);
            sObj.setName(name);            
        end
        function setUnits(sObj, xUnits, yUnits)
            %setUnits(sObj, xUnits, yUnits)
            %Set the axis units.
            %Same as calling sObj.setXUnits(xUnits) and
            %sObj.setYUnits(yUnits) separately.
            %yUnits is optional argument. If it is not specified, the this
            %function behaves like setXUnits.
            if(nargin==3)
                if(isa(yUnits,'char'))
                    sObj.setYunits(yUnits);
                end
            end
            if(nargin>=2)
                if(isa(xUnits,'char'))
                    sObj.setXunits(xUnits);
                end
            end
        end
        function setXUnits(sObj, units)
            %setXUnits(sObj, units)
            %Sets the units of the x-axis
            if(isa(units, 'char'))
                sObj.xunits = units;
            end
        end
        function setYUnits(sObj, units)
            %setYUnits(sObj, units)
            %Sets the units of the y-axis
            if(isa(units,'char'))
                sObj.yunits = units;
            end
        end
        function setSampleRate(sObj, sampleRate)
            % setSampleRate(sObj, sampleRate)
            % sets the current sampleRate of the object to rate specified
            if(sObj.sampleRate~=sampleRate)
                if(~(floor(sampleRate*200)/200==floor(sObj.sampleRate*200)/200)) %Compare to 2 decimal places (finite precision has caused errors 500.000001 ~= 500.00000x
                    if(sampleRate>sObj.sampleRate)
                        %fprintf(strcat('SignalObj,',sObj.name',', upsampled to:',num2str(sampleRate)));
                    else
                        %fprintf(strcat('SignalObj,',sObj.name',', downsampled to:',num2str(sampleRate)));
                    end
                    sObj.resampleMe(sampleRate);
                end
            end
        end
        function setDataLabels(sObj,dataLabels)
            %setDataLabels(sObj,dataLabels)
            %sets the labels for each of the components of the SignalObj.
            %if sObj has only a single component, then dataLabels can be a
            %string. Otherwise, dataLabels must be a cell with the same
            %dimensions as sObj. dataLabels{i} specifies the string for the
            %ith component of sObj.
            if(~isempty(dataLabels))
                if(isa(dataLabels,'char'))
                    if(sObj.dimension==1)
                        sObj.dataLabels{1}=dataLabels;
                    else
                        display('Adding single dataLabel to a SignalObj with more that 1 dimension. All dimensions have same label now!');
                        for i=1:sObj.dimension
                            sObj.dataLabels{i} = dataLabels;
                        end
                    end

                elseif(isa(dataLabels,'cell'))
                    if(length(dataLabels)==sObj.dimension)
                        %ind=sObj.findIndFromDataMask;
                        %for i=ind
                        %    sObj.dataLabels{i} = dataLabels{i};
                        %end     
                        sObj.dataLabels = dataLabels;
                    else
                        error('Need the number of labels to match the number of dimensions of the SignalObj');
                    end
                end
            end
        end
        function setMinTime(sObj,minTime,holdVals)
            %setMinTime(sObj,minTime,holdVals)
            %sets the minimun value of the time vector to minTime. If
            %minTime>min(sObj.time) then the data before minTime will be
            %ignored. If minTime < min(sObj.time) then the time vector is
            %extended at the current sampleRate to minTime. 
            %holdVals: 1 or 0. If not specifed, defaults to 0. If
            %holdVals=1, then the value at min(sObj.time) is extended to
            %the new minTime. Otherwise, the added time is padded with
            %zeros.
            if(nargin<3)
                holdVals=0;
            end
            
            if(nargin<2)
                minTime=sObj.time(1);
            end
            timeVec=sObj.getTime;
            if(minTime<min(timeVec))
                maxTime=max(timeVec);
                newTime=minTime:1/sObj.sampleRate:maxTime;
                newTime=newTime';
                numSamples = length(newTime)-length(timeVec);
                if(holdVals==1)
                    newData=[ones(numSamples,1)*sObj.data(1,:);sObj.data];
                else
                    newData=[zeros(numSamples,sObj.dimension);sObj.data];
                end
                sObj.data=newData;
                sObj.time=newTime;
                sObj.minTime=min(sObj.time);
            elseif(min(timeVec)==minTime)
                    %do nothing
            else
                startIndex = sObj.findNearestTimeIndex(minTime);
                sObj.time=sObj.time(startIndex:end);
                sObj.data=sObj.data(startIndex:end,:);
            end
            sObj.minTime=min(sObj.time);
        end
        function setMaxTime(sObj,maxTime, holdVals)
            %setMaxTime(sObj,maxTime,holdVals)
            %sets the maximum value of the time vector to maxTime. If
            %maxTime<max(sObj.time) then the data after maxTime will be
            %ignored. If maxTime > max(sObj.time) then the time vector is
            %extended at the current sampleRate to maxTime. 
            %holdVals: 1 or 0. If not specifed, defaults to 0. If
            %holdVals=1, then the value at min(sObj.time) is extended to
            %the new minTime. Otherwise, the added time is padded with
            %zeros.
            if(nargin<3)
                holdVals=0;
            end
            
            if(nargin<2)
                maxTime=sObj.time(end);
            end
            timeVec=sObj.getTime;
            if(max(timeVec)<maxTime)
                minTime=min(timeVec);
                newTime=linspace(minTime,maxTime,(sObj.sampleRate)*(maxTime-minTime)+1);
                newTime = newTime';
                numSamples = length(newTime)-length(timeVec);
                if(holdVals==1)
                    newData=[sObj.data;ones(numSamples,1)*sObj.data(end,:)];
                else
                    newData=[sObj.data;zeros(numSamples,sObj.dimension)];
                end
                
                sObj.data=newData;
                sObj.time=newTime;
                sObj.maxTime=max(sObj.time);
            elseif(max(timeVec)==maxTime)
                    %do nothing
                
            else
                endIndex = sObj.findNearestTimeIndex(maxTime);
                sObj.time=sObj.time(1:endIndex);
                sObj.data=sObj.data(1:endIndex,:);
            end
            sObj.maxTime=max(sObj.time);
        end
        function setPlotProps(sObj, plotProps,index)
            %setPlotProps(sObj, plotProps,index)
            %if index is not specified:
            %   - plotProps is a cell with sObj.dimension elements, then plotProps{i} specifies a string
            %     that will be used to plot the ith component of sObj. 
            %   - plotProps is a string, then the string will be used to
            %     plot all of the components of sObj.
            %
            %if index is specified and index is within range of the number
            %of components of the signal:
            %   - plotProps is a cell of length 1 then the property is
            %     applied to the component specified by the index
            %   - plotProps is a string, the property is applied to the
            %     component specified by the index.
            if(nargin<=2)
                if(isa(plotProps,'cell'))
                    if(length(plotProps) == sObj.dimension)
                        for i=1:sObj.dimension
                            sObj.plotProps{i} = cell2str(plotProps{i});
                        end
                    elseif(length(plotProps)==1)
                        for i=1:sObj.dimension
                            sObj.plotProps{i} = cell2str(plotProps);
                            display('Index not specified. All dimensions set to have same plotting properties');
                        end
                    else
                        error('Index not specified and more than 1 plotProp specified. Need to number of plotProps same as sObj.dimension or length 1');
                    end
                elseif(isa(plotProps,'char'))
                     for i=1:sObj.dimension
                        sObj.plotProps{i} = cell2str(plotProps);
                     end
                     display('All dimensions set to have same plotting properties')
                end
               
            else
                if(isa(plotProps,'cell') && length(plotProps)==1)
                    if(index>0 && index<=sObj.dimension)
                        sObj.plotProps{index} = plotProps{:};
                    else
                        error('Index out of bounds during setPlotProps');
                    end
                elseif(isa(plotProps,'char'))
                     if(index>0 && index<=sObj.dimension)
                        sObj.plotProps{index} = plotProps;
                    else
                        error('Index out of bounds during setPlotProps');
                    end
                    
                end
            end
        end  
        function setMask(sObj, mask)
            %setMask(sObj, mask)
            % if called with no arguments, all the components of the signal
            % are masked. No data will be visible.
            % mask: either a set of indices or a cell array of characters
            % indicating which signal components are to remain visible.
            if(nargin<2)
                mask=zeros(1,sObj.dimension);
                sObj.setDataMask(mask);
                return;
            end
            %mask is either a set of indices or names;
            if(isa(mask,'cell'))
                if(isa(mask{1},'char'))
                    sObj.setMaskByLabels(mask);
                else
                    error('Mask cells must contains strings!');
                end
            elseif(isa(mask,'double'))
                sObj.setMaskByInd(mask);
            else
                error('Can only set datamask with strings or indices')
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function tVec       = getTime(sObj)
            % tVec       = getTime(sObj)
            % returns the time vector of the Signal Obj
            tVec=sObj.time;
        end
        function data       = getData(sObj)
            % data       = getData(sObj)
            % Returns the signal data as a matrix. If masks are set, then
            % only the components that are visible will be returned. Each
            % column corresponds to each component of the SignalObj that is
            % visible. The columns are in the same order as the dataLabels.
            data=sObj.dataToMatrix;
        end
        function [t,d]      = getOriginalData(sObj)
            % [t,d]      = getOriginalData(sObj)
            % SignalObjs have memory. The original data and time vectors
            % are stored even when the signal is resamples, windowed, etc.
            % This commands returns the original data used to create the
            % SignalObj
            t=sObj.originalTime;
            d=sObj.originalData;
        end
        function s          = getOrigDataSig(sObj)
            % s          = getOrigDataSig(sObj)
            % same as getOriginalData, except that a SignalObj containing the original data is returned.
            
            [time,data]=sObj.getOriginalData;
            name=sObj.name;
            xlabelval=sObj.xlabelval; 
            xunits=sObj.xunits;
            yunits=sObj.yunits; 
            dataLabels=sObj.dataLabels;
            plotProps=sObj.plotProps;
            evalstring = strcat('s=',class(sObj),'(time, data,name, xlabelval, xunits, yunits,dataLabels,plotProps)');
            eval(evalstring);
            %s = SignalObj(time, data,name, xlabelval, xunits, yunits,dataLabels,plotProps);
        end
        function val        = getValueAt(sObj,x)
            %val        = getValueAt(sObj,x)
            %returns a row vector of length sObj.dimension corresponding to
            %the values of the signal evaluated at time=x
            
            %ind=sObj.findNearestTimeIndices(x);
            %val=sObj.data(ind,:);
            
            [l,w]=size(x);
            if(w>l)
                x=x';
            end 
           
%             val = interp1(sObj.time,sObj.data,x,'spline',0); %extrapolate to zero
            val = interp1(sObj.time,sObj.data,x,'nearest',0); %extrapolate to zero
%             if(any(isnan(sObj.data)))
%                 pause
%             end
%             if(sObj.dimension==1)
%                 val=val';
%             end
           
        end
        function PropsStr   = getPlotProps(sObj,index)
            %PropsStr   = getPlotProps(sObj,index)
            %Returns the string correspond to the plotting properties of
            %the SignalObj component corresponding to index
            if(index>0 && index<=sObj.dimension)
                PropsStr = cell2str(sObj.plotProps{index});
            else
                error('index is out of bounds!');
            end
        end
        function indices    = getIndicesFromLabels(sObj,label)
            %indices    = getIndicesFromLabels(sObj,label)
            %Returns a cell array if the label appears for various point
            %in the SignalObj. indices{i} contains all the the indices
            %corresponding to label{i} if label is a cell-array or label if
            %it is a string.
            %Returns an array if the SignalObj label appears only once in the
            %SignalObj
            if(isa(label,'cell'))
                indices =cell(1,length(label));
                numInd  =zeros(1,length(label));
                for i=1:length(label)
                    tempInd = sObj.getIndexFromLabel(label{i});
                    if(~isempty(tempInd))
                        numInd(i) = length(tempInd);
                        indices{i}=tempInd;
                    else
                        error('Label does not exist!');
                    end
                     
                end
            elseif(isa(label,'char'))
                indices = sObj.getIndexFromLabel(label);
                numInd(1) = length(indices);
            end
            

            if(max(numInd)==1)      %For backwards compatibility if assuming only on index per label
                if(isa(indices,'cell'))
                    for i=1:length(numInd)

                        tempInd(i) = indices{i};
                    end
                    indices = tempInd;
                end
            end
        end  
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Operand Definitions and other mathematical operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function s3 = plus(s1,s2)
            % s3 = plus(s1,s2)
            % Adds two signals
            if(isa(s1,'SignalObj') && isa(s2,'SignalObj'))
                % What if s2 is a constant or double?
                if(s1.dimension == s2.dimension)
                    [s1c,s2c] = makeCompatible(s1,s2);
                    s3=s1c.copySignal;
                    s3.data = s1c.data+s2c.data;
                    
                    for i=1:length(s3.dataLabels)
                        if(~s2c.areDataLabelsEmpty && ~isempty(s2c.dataLabels{i}))
                            if(strcmp(s2c.dataLabels{i}(1),'-'))
                                s3.dataLabels{i} = [s1c.dataLabels{i} '-' s2c.dataLabels{i}(2:end)];
                            else
                                s3.dataLabels{i} = [s1c.dataLabels{i} '+' s2c.dataLabels{i}];
                            end
                        end
                    end
                else
                    error('Can only add signals if they have the same dimension');
                end
            elseif(isa(s1,'double') || isa(s2,'double'))
                if(isa(s1,'double'))
                  s3=s2.copySignal;
                  [l,w] = size(s1);
                  if(w==s3.dimension && l==1)
                      s3.data = s3.data+ones(length(s3.data),1)*s1;
                      for i=1:length(s3.dataLabels)
                            if(sign(s1(i))==-1)
                                s3.dataLabels{i} = [s2.dataLabels{i} '-' num2str(abs(s1(i)))];
                            else
                                s3.dataLabels{i} = [s2.dataLabels{i} '+' num2str(abs(s1(i)))];
                            end
                      end
                  else
                      s3.data = s3.data+s1; 
                      %dont modify dataLabels since s1 is a matrix;
                             
                      %for i=1:length(s3.dataLabels)
                      %      s3.dataLabels{i} = [s2.dataLabels{i} '+' num2str(s1(i))];
                      %end
                  end
                else
                  s3=s1.copySignal;
                  [l,w] = size(s2);
                  %size(s3.data)
                  if(w==s3.dimension && l==1)
                      s3.data = s3.data+ones(length(s3.data),1)*s2;
                      for i=1:length(s3.dataLabels)
                            if(sign(s2(i))==-1)
                                s3.dataLabels{i} = [s1.dataLabels{i} '-' num2str(abs(s2(i)))];
                            else
                                s3.dataLabels{i} = [s1.dataLabels{i} '+' num2str(abs(s2(i)))];
                            end
                      end
                  else
                      
                      s3.data = s3.data+s2; 
                      %dont modify dataLabels since s2 is a matrix;
                      %for i=1:length(s3.dataLabels)
                      %      s3.dataLabels{i} = [s1.dataLabels{i} '+' num2str(s2)];
                      %end
                  end
                end
            else
                error('only Signals or doubles are currently supported');
            end
        end        
        function s3 = minus(s1,s2)
            % s3 = minus(s1,s2)
            % Subtracts two signals
            s3=plus(s1,-s2);
        end        
        function s3 = uplus(s1) %+s1
            % s3 = uplus(s1)
            % Multiplies signal by +1
            s3=s1.copySignal;
        end        
        function s3 = uminus(s1) %-s1
            % s3 = uminus(s1) %-s1
            % Multiplies the signal by -1.
            % dataLabels are updated to reflect this.
            s3=s1.copySignal;
            s3.data=-s3.data;
            for i=1:length(s3.dataLabels)
                s3.dataLabels{i} = strcat('-',s1.dataLabels{i});
            end
        end 
        function s3 = power(s1,exponent)
           if(isa(exponent,'double'))
               s3=s1.copySignal;
               s3.data=s3.data.^exponent;
           else
               error('Exponent should be a double');
           end
        end
        
        function s3 = sqrt(s1)
           s3=s1.copySignal;
           s3.data=sqrt(s3.data);
            
        end
        function s3 = times(s1,s2) %s1.*s2
            %s3 = times(s1,s2)
            %Multiplies each sample in s1 with each sample of s2.
            %s1 or s2 can be doubles
            if(isa(s1,'SignalObj') && isa(s2,'SignalObj'))
               if(s1.dimension == s2.dimension)
                   [s1c,s2c] = makeCompatible(s1,s2);
                   s3 = s1c;
                   s3.data = s1c.data.*s2c.data;
                   %can multiply units 
               else
                   error('can only multiply signals with same dimension');
               end
            elseif(isa(s1,'double') || isa(s2,'double'))
                if(isa(s1,'double'))
                  s3=s2.copySignal;
                  [l,w] = size(s1);
                  if(w==s3.dimension && l==1)
                      s3.data = s3.data.*(ones(length(s3.data),1)*s1);
                  else
                      s3.data = s3.data.*s1;
                  end
                else
                  s3=s1.copySignal;
                  [l,w] = size(s2);
                  if(w==s3.dimension && l==1)
                      s3.data = s3.data.*(ones(length(s3.data),1)*s2);
                  else
                      s3.data = s3.data.*s2;
                  end
                end 
            end    
        end        
        function s3 = mtimes(s1,s2) %s1*s2
            %Matrix multiplication of two signals
            %Needs work
            %If s1 and s2 are signals, same as times(s1,s2)
            if(isa(s1,'SignalObj') && isa(s2,'SignalObj'))
                   %[s1c,s2c] = makeCompatible(s1,s2);
                   s3 = s1.copySignal;
                   s3.data = s1.data.*s2.data;
                   %can multiply units 
           elseif(isa(s1,'double') || isa(s2,'double'))
                if(isa(s1,'double'))
                  s3=s2.copySignal;
                  s3.data = (s1*s3.data')';
                else
                  s3=s1.copySignal;
                  s3.data = (s3.data'*s2)';
                end 
            end
        end        
        function s3 = rdivide(s1,s2)
            if(isa(s1,'SignalObj') && isa(s2,'SignalObj'))
               if(s1.dimension == s2.dimension)
                   [s1c,s2c] = makeCompatible(s1,s2);
                   s3 = s1c;
                   s3.data = s1c.data./s2c.data;
                   %can multiply units 
               else
                   error('can only multiply signals with same dimension');
               end
            elseif(isa(s1,'double') || isa(s2,'double'))
                if(isa(s1,'double'))
                  s3=s2.copySignal;
                  s3.data = s1./s3.data;
                else
                  s3=s1.copySignal;
                  s3.data = s3.data./s2;
                end 
            end
        end        
        function s3 = ldivide(s1,s2)
            if(isa(s1,'SignalObj') && isa(s2,'SignalObj'))
               if(s1.dimension == s2.dimension)
                   [s1c,s2c] = makeCompatible(s1,s2);
                   s3 = s1c;
                   s3.data = s1c.data.\s2c.data;
                   %can multiply units 
               else
                   error('can only multiply signals with same dimension');
               end
            elseif(isa(s1,'double') || isa(s2,'double'))
                if(isa(s1,'double'))
                  s3=s2.copySignal;
                  s3.data = s1.\s3.data;
                else
                  s3=s1.copySignal;
                  s3.data = s3.data.\s2;
                end 
            end
        end          
        function s3 = ctranspose(s1)
            s3=s1.copySignal;
            s3.data=s3.data.';
            [l,w]=size(s3.data);
            s3.dimension=w;
        end        
        function s3 = transpose(s1)
            s3=s1.copySignal;
            s3.data=s3.data.';
            [l,w]=size(s3.data);
            s3.dimension=w;
        end
        function s3 = derivative(sObj)
            %Computes derivative of each component of the SignalObj with
            %respect to the x-axis variable.
%             B=[1 -1]*sObj.sampleRate;
%             A=1;
%             s3=sObj.filter(B,A);
            s3=sObj.copySignal;
            tData=diff(s3.data)*s3.sampleRate;
            tData=[zeros(1,s3.dimension); tData];
            s3.data=tData;
            
            s3.setYUnits(strcat('\frac{',s3.yunits,'}{',s3.xunits,'}'));
            denomstr = strcat('d', s3.xlabelval(1));
            
            s3.setName(strcat('\frac{d}{',denomstr,'}',s3.name));
            for i=1:s3.dimension
                if(~strcmp(sObj.dataLabels{i},''))
                    s3.dataLabels{i}= strcat('\frac{d}{',denomstr,'}',s3.dataLabels{i});
                end
            end
        end        
        function val = derivativeAt(sObj,x0)
            % val = derivativeAt(sObj,x0)
            %computes the derivative of the Signal at x0. Returns a row
            %vector of length equal to sObj.dimension.
            sTemp = sObj.derivative;
            val = sTemp.getValueAt(x0);
        end        
        function s3 = integral(sObj,t0,tf)
            %s3 = integral(sObj,t0,tf)
            %computes the integral of the signal in the window from t0 to tf.
            %if tf is not specified, sObj.maxTime is used.
            %if t0 is not specified, sObj.minTime is used.
            % Both t0 and tf are optional but t0 must be specified if tf is
            % to be specified (e.g. cant specified tf only)
            % Data labels are updated with latex notation for integral.
            % the value of the returned signal at time t is the value of
            % the integral from minTime to t.
            if(nargin<3)
                tf=sObj.maxTime;
            end
            if(nargin<2)
                t0=sObj.minTime;
            end
            
            %y[n] = y[n-1] + x[n]*deltaT
            B=1*1/sObj.sampleRate;
            A=[1 -1];
            s3=sObj.getSigInTimeWindow(t0,tf);
            s3=s3.filter(B,A);
            s3.setYUnits(strcat(s3.yunits,'*',s3.xunits));
            dtstr = strcat(' d','\tau');
            s3.setName(['\int_',num2str(s3.minTime),'^',s3.xlabelval(1),'\!\!{',[s3.name dtstr],'}']);
            if(~sObj.areDataLabelsEmpty)
                for i=1:s3.dimension
                    if(~strcmp(sObj.dataLabels{i},''))
                        s3.dataLabels{i}= ['\int_',num2str(s3.minTime),'^',s3.xlabelval(1),'\!\!{',[s3.dataLabels{i} dtstr],'}'];
                    else
                        s3.dataLabels{i} = '';
                    end
                end
            end
        end          
        function s3 = filter(sObj, B,A)
            %s3 = filter(sObj, B,A)
            %applies the discrete filter specified by B and A to the sObj
            %data. Same as running filter(B,A,sObj.dataToMatrix). dataMasks
            %are ignores so that the signal dimensionality does not change.
            s3=sObj.copySignal;
            s3.data = filter(B,A,s3.data);
        end        
        function s3 = filtfilt(sObj,B,A)
            %s3 = filtfilt(sObj,B,A) 
            %applies the discrete filter specified by B and A to the sObj
            %data using filtfilt. Same as running filtfilt(B,A,sObj.dataToMatrix). dataMasks
            %are ignores so that the signal dimensionality does not change.
            s3=sObj.copySignal;
            s3.data = filtfilt(B,A,s3.data);
        end
        function [s1c,s2c] = makeCompatible(s1,s2,holdVals)
            %[s1c,s2c] = makeCompatible(s1,s2,holdVals)
            %returns two signals that copies of the original signals but
            %that have been resampled or extended in time so that the time
            %axis of both signals in the same. This is done before most
            %mathmatical operations to make sure that the signals have the
            %same support.
            % holdVals = 1 makes the signals keep their endpoint values if
            % they are extended in time. holdVals is an optional argument.
            if(nargin<3)
                holdVals=0;
            end
            if(s1.minTime~=s2.minTime || s1.maxTime~=s2.maxTime || s1.sampleRate ~=s2.sampleRate)
                s1c = s1.copySignal; s2c = s2.copySignal;
                minTime=min(s1c.minTime,s2c.minTime);
                maxTime=max(s1c.maxTime,s2c.maxTime);
                sampleRate=max(s1c.sampleRate,s2c.sampleRate);
                s1c.setSampleRate(sampleRate); s2c.setSampleRate(sampleRate);
                s1c.setMinTime(minTime,holdVals);       s2c.setMinTime(minTime,holdVals);
                s1c.setMaxTime(maxTime,holdVals);       s2c.setMaxTime(maxTime,holdVals);

                %pause
    %             for i=1:s2c.dimension
                 %if(max(s2c.time-s1c.time)>0)
                    data = interp1(s2c.time,s2c.data,s1c.time,'nearest',0);
                 %else
                 %   data = s2c.data;
                 %end
    %             end
                 s2c.time = s1c.time;
                 [nrows,ncolumns] = size(data);
                 if(nrows>ncolumns)
                    s2c.data = data;
                 else
                    s2c.data = data';
                 end
            else
                s1c = s1;
                s2c = s2;
            end
                
        end
        function s = abs(sObj)
            absData=abs(sObj.data);
            [nrows,ncolumns]=size(absData);
            name =  ['|', sObj.name '|'];
            plotProps  = sObj.plotProps;
            if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
%                     plotProps  = sObj.plotProps;
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('|',sObj.dataLabels{i},'|');
                    end
                evalstring = strcat('s=',class(sObj),'(sObj.time, absData,name,sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels,plotProps);');
            else
                
                evalstring = strcat('s=',class(sObj),'(sObj.time, absData,name,sObj.xlabelval, sObj.xunits,sObj.yunits,[],plotProps);');
            end
            eval(evalstring);            
        end
        
        function s = log(sObj)
            logData=log(sObj.data);
            [nrows,ncolumns]=size(logData);
            name =  ['ln(' sObj.name ')'];
            yunits = ['ln(' sObj.yunits ')'];
            plotProps  = sObj.plotProps;
            if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
%                     plotProps  = sObj.plotProps;
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('ln(',sObj.dataLabels{i},')');
                    end
                evalstring = strcat('s=',class(sObj),'(sObj.time, logData,name,sObj.xlabelval, sObj.xunits,yunits,dataLabels,plotProps);');
            else
                
                evalstring = strcat('s=',class(sObj),'(sObj.time, logData,name,sObj.xlabelval, sObj.xunits,yunits,[],plotProps);');
            end
            eval(evalstring);            
        end
        
        function m=median(sObj,varargin)
            %m=median(sObj,varargin)
            %Computes the column-wise median of SignalObj data. Returns a
            %signal with the corresponding median values
            %same as calling median(sObj.dataToMatrix,varargin)
            %Additional parameters are passed to the matlab median function.
            %Default computes median of each signal component across time.
            %mean(sObj,2) computes median value of the SignalObj at each
            %point in time.
            mdata=median(sObj.data,varargin{:});
            [nrows,ncolumns]=size(mdata);
            if( (nrows==length(sObj.time)) && (ncolumns==1) ) %mean across dimensions
                name =  ['median(', sObj.name ')'];
                evalstring = strcat('m=',class(sObj),'(sObj.time, mdata,name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                eval(evalstring);
            elseif( (nrows==1) && (ncolumns == sObj.dimension) )  %mean of each dimension
                if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('median(',sObj.dataLabels{i},')');
                    end
                     name =  ['median(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata], ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);
                else
                     name =  ['median(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata],  ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits);
                end
                
            
            end
           
            
            
        end
        function m=mode(sObj,varargin)
            %m=mode(sObj,varargin)
            %Computes the column-wise mode of SignalObj data. Returns a
            %signal with the corresponding mode values
            %same as calling mode(sObj.dataToMatrix,varargin)
            %Additional parameters are passed to the matlab mode function.
            %Default computes mode of each signal component across time.
            %mean(sObj,2) computes mode value of the SignalObj at each
            %point in time.
            mdata=mode(sObj.data,varargin{:});
            [nrows,ncolumns]=size(mdata);
            if( (nrows==length(sObj.time)) && (ncolumns==1) ) %mean across dimensions
                name =  ['mode(', sObj.name ')'];
                evalstring = strcat('m=',class(sObj),'(sObj.time, mdata,name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                eval(evalstring);
            elseif( (nrows==1) && (ncolumns == sObj.dimension) )  %mean of each dimension
                if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('mode(',sObj.dataLabels{i},')');
                    end
                     name =  ['mode(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata], ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);
                else
                     name =  ['mode(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata],  ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits);
                end
            end
        end
        
        
        function m=mean(sObj,varargin)
            %m=mean(sObj,varargin)
            %Computes the column-wise mean of SignalObj data. Returns a
            %signal with the corresponding mean values
            %same as calling mean(sObj.dataToMatrix,varargin)
            %Additional parameters are passed to the matlab mean function.
            %Default computes mean of each signal component across time.
            %mean(sObj,2) computes mean value of the SignalObj at each
            %point in time.
            mdata=mean(sObj.data,varargin{:});
            [nrows,ncolumns]=size(mdata);
            if( (nrows==length(sObj.time)) && (ncolumns==1) ) %mean across dimensions
                name =  ['\mu(', sObj.name ')'];
                evalstring = strcat('m=',class(sObj),'(sObj.time, mdata,name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                eval(evalstring);
            elseif( (nrows==1) && (ncolumns == sObj.dimension) )  %mean of each dimension
                if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('\mu(',sObj.dataLabels{i},')');
                    end
                     name =  ['\mu(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata], ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);
                else
                     name =  ['\mu(', sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [mdata;mdata],name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [mdata;mdata],  ['Mean of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits);
                end
            end
        end
        function m=std(sObj,varargin)
            %m=std(sObj,varargin)
            %Computes the column-wise std of SignalObj data. Returns a
            %signal with the corresponding mean values
            %same as calling std(sObj.dataToMatrix,varargin)
            %Additional parameters are passed to the matlab std function.
            %Default computes std of each signal component across time.
            %std(sObj,[],2) computes std value of the SignalObj at each
            %point in time.
            stdData=std(sObj.data,varargin{:});
            [nrows,ncolumns]=size(stdData);
            if( (nrows==length(sObj.time)) && (ncolumns==1) ) %mean across dimensions
                     name=['\sigma(' sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'(sObj.time, stdData,name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                     eval(evalstring);
                     %m=SignalObj(sObj.time, stdData, name,sObj.xlabelval, sObj.xunits,sObj.yunits); 
            elseif( (nrows==1) && (ncolumns == sObj.dimension) )  %mean of each dimension
                if(~sObj.areDataLabelsEmpty)
                    dataLabels = cell(size(sObj.dataLabels));
                    for i=1:sObj.dimension
                        dataLabels{i} = strcat('\sigma(',sObj.dataLabels{i},')');
                    end
                     name=['\sigma(' sObj.name ')'];
                     evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [stdData;stdData],name,sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);');
                     eval(evalstring);
                     %m=SignalObj([sObj.time(1); sObj.time(end)], [stdData;stdData], ['Std. Dev. of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits,dataLabels);
                else
                    name=['\sigma(' sObj.name ')'];
                    evalstring = strcat('m=',class(sObj),'([sObj.time(1); sObj.time(end)], [stdData;stdData],name,sObj.xlabelval, sObj.xunits,sObj.yunits);');
                    eval(evalstring);
                    %m=SignalObj([sObj.time(1); sObj.time(end)], [stdData;stdData],  ['Std. Dev. of ' sObj.name],sObj.xlabelval, sObj.xunits,sObj.yunits);
                end
            end
            
        end
        function [m, index,time]=max(sObj,varargin)
            %[m,time]=max(sObj,varargin)
            %calls matlab max function with the SignalObj's data.
            %Additional parameters are passed to the matlab function
            [m,index]=max(sObj.data,varargin{:});
            time = sObj.time(index);
            
        end
        function [m,index,time]=min(sObj,varargin)
            %m=min(sObj,varargin)
            %calls matlab min function with the SignalObj's data.
            %Additional parameters are passed to the matlab function
            [m,index]=min(sObj.data,varargin{:});
            time = sObj.time(index);
        end
        function periodogram = periodogram(sObj)
            %periodogram = periodogram(sObj)
            % computes the periodogram of each component of the SignalObj.
            % Calls matlab periodogram with each component. 
            %   fs = sObj.sampleRate;
            %   NFFT = 1024
            fs = sObj.sampleRate;  % Sampling frequency
            periodogram = cell(1,sObj.dimension);
            for i=1:sObj.dimension
                xn = sObj.data(:,i);
                Hs = spectrum.periodogram('rectangular');
                switch sObj.dimension
                    case 2
                        subplot(1,2,i)
                    case 3
                        subplot(1,3,i)
                    case 4
                        subplot(2,2,i)
                    case 5
                        subplot(3,2,i)
                    case 6
                        subplot(3,2,i)
                    otherwise
                         h=gcf;figure(h);
                end
                periodogram{i}=psd(Hs,xn,'Fs',fs,'NFFT',1024);
                h=periodogram{i}.plot;legend(h, sObj.dataLabels{i});
            end
        end
        function mtmSpec = MTMspectrum(sObj,NW,NFFT,Pval)
            %mtmSpec = MTMspectrum(sObj,NW,NFFT,Pval)
            %computes Multi-taper spectral estimate of each component of
            %the SignalObj with defaults:
            %NW=4, NFFT=[], Pval=.95
            % Defaults can be changed by passing in additional arguments
            if(nargin<4)
                Pval=.95;
            end
            if(nargin<3)
                NFFT=[];
            end
            if(nargin<2)
                NW=4;
            end
            
            Fs=sObj.sampleRate;
            mtmSpec = cell(1,sObj.dimension);
            for i=1:sObj.dimension
                xn=sObj.data(:,i);
                [Pxx,Pxxc,f] = pmtm(xn,NW,NFFT,Fs,Pval);
                hpsd = dspdata.psd([Pxx Pxxc],'Fs',Fs);
                mtmSpec{i} = hpsd;
                switch sObj.dimension
                    case 2
                        subplot(1,2,i)
                    case 3
                        subplot(1,3,i)
                    case 4
                        subplot(2,2,i)
                    case 5
                        subplot(3,2,i)
                    case 6
                        subplot(3,2,i)
                    otherwise
                         h=gcf;figure(h);
                end
                str1=strcat(num2str(Pval*100), '% Conf. Int.');
                h=plot(hpsd); legend(h, sObj.dataLabels{i},strcat('-',str1),strcat('+',str1));
            end
        end
        function h = spectrogram(sObj)
            t=sObj.time;             % 2 secs @ 1kHz sample rate
            x=sObj.data;             % Start @ DC, cross 150Hz at t=1sec 
            F = 0:.1:100;
            
            clear y f t p;
            for i=1:sObj.dimension;
                figure;
                [y{i},f{i},t{i},p{i}] = spectrogram(x(:,i),256,250,F,sObj.sampleRate,'yaxis'); 
                surf(t{i},f{i},10*log10(abs(p{i})),'EdgeColor','none');   
                axis xy; axis tight; colormap(jet); view(0,90);
                xlabel('Time');
                ylabel('Frequency (Hz)');
            end
            
        end
        
        
        function sxCorr= xcorr(s1,s2,varargin)
            if(nargin<2)
                s2=s1;
            end
            
             [s1c,s2c] = makeCompatible(s1,s2);
             if(~isempty(varargin))
                [tempC tempLags] =xcorr([s1c.data, s2c.data],varargin{1});
             else
                 [tempC tempLags] =xcorr([s1c.data, s2c.data]);
             end
             index=[];
             for i=1:s1c.dimension
                 offset(i)=(i-1)*(s1c.dimension+s2c.dimension);
                 index = [index, (offset(i)+s1c.dimension+1):(offset(i)+s1c.dimension+s2c.dimension)];
             end
             sum=0;
             dataLabels = cell(1,length(index));
             for i=1:s1c.dimension
                 for j=1:s2c.dimension
                     sum = sum+1;
                     dataLabels{sum} = strcat('corr(',s1c.dataLabels{i},',',s2c.dataLabels{j}, ']');
                 end
             end
             
             M=length(s1c.data);
             if(nargin<2)
                data=tempC(M-1:end,index);
                lags=tempLags(M-1:end)./s1c.sampleRate;
             else
                 data=tempC(1:end,index);
                lags=tempLags(1:end)./s1c.sampleRate;
             end
             name = [ 'corr(' s1.name ',' s2.name ')'];
             evalstring=strcat('sxCorr=',class(s1),'(lags,data,name,''\Delta \tau'',''s'',dataLabels);');
             eval(evalstring);    
        end        
        function sxCov = xcov(s1,s2,varargin)
             if(nargin<2)
                s2=s1;
             end
             [s1c,s2c] = makeCompatible(s1,s2);
             if(~isempty(varargin))
                [tempC tempLags] =xcov([s1c.data, s2c.data],varargin{1});
             else
                 [tempC tempLags] =xcov([s1c.data, s2c.data]);
             end
             index=[];
             for i=1:s1c.dimension
                 offset(i)=(i-1)*(s1c.dimension+s2c.dimension);
                 index = [index, (offset(i)+s1c.dimension+1):(offset(i)+s1c.dimension+s2c.dimension)];
             end
             sum=0;
             dataLabels = cell(1,length(index));
             for i=1:s1c.dimension
                 for j=1:s2c.dimension
                     sum = sum+1;
                     dataLabels{sum} = strcat('cov(',s1c.dataLabels{i},',',s2c.dataLabels{j}, ']');
                 end
             end
             
             M=length(s1c.data);
             if(nargin<2)
                data=tempC(M-1:end,index);
                lags=tempLags(M-1:end)./s1c.sampleRate;
             else
                 data=tempC(1:end,index);
                 lags=tempLags(1:end)./s1c.sampleRate;
             end
             name = [ 'cov(' s1.name ',' s2.name ')'];
             evalstring=strcat('sxCov=',class(s1),'(lags,data,name,''\Delta \tau'',''s'',dataLabels);');
             eval(evalstring);    
             %sxCov=SignalObj(lags,data,name,'lags','n',dataLabels);

                 
         end
        
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function mergedSig = merge(sObj,varargin)
            %mergedSig = merge(sObj,varargin)
            %merges the data in any signals passed in with the data in sObj.
            %if one of the input arguments is a double then it is assumed
            %to be the value of holdVals to be used when making the signals
            %compatible. See makeCompatible(s1,s2)
            numToMerge=0;
            holdVals=0;
            for i=1:length(varargin)
                if(isa(varargin{i},'SignalObj'))
                    numToMerge = numToMerge+1;
                end
                if(isa(varargin{i},'double')&&i==length(varargin)) %expect only one double at the end
                    holdVals = varargin{i};
                end
            end

            if(numToMerge == 1)
                
                [s1c, s2c] = sObj.makeCompatible(varargin{1},holdVals);
                 data=[s1c.dataToMatrix, s2c.dataToMatrix];
                 dataLabels = cell(s1c.dimension+s2c.dimension,1);
                 for i=1:s1c.dimension
                        dataLabels{i} = s1c.dataLabels{i}; 
                 end
                 for i=1:s2c.dimension
                        dataLabels{s1c.dimension+i}=s2c.dataLabels{i};
                 end
                 name = s1c.name;%, ',', s2c.name]; %name of the original is kept
                 
                 evalstring = strcat('mergedSig = ',class(sObj),'(s1c.time, data, name, s1c.xlabelval,s1c.xunits, s1c.yunits, dataLabels);');
                 eval(evalstring);
            else
                  mergedSig = sObj.merge(varargin{1},holdVals);
                  for i=2:numToMerge
                        mergedSig = mergedSig.merge(varargin{i},holdVals);
                  end
            end
        end        
        function sigOut=copySignal(sigIn)
            %sigOut=copySignal(sigIn)
            %sigOut is a copy of the original signal sigIn.
            %This function is necessary because SignalObj's are handles and
            %as such, setting s2=s1, and then changing s1, causes changes
            %in s2 as well. To avoid this, make a copy of s1 and set it as
            %s2. (eg. s2=s1.copySignal; )
                time=sigIn.time; 
                data=sigIn.data;
                name=sigIn.name;
                xlabelval=sigIn.xlabelval; 
                xunits=sigIn.xunits;
                yunits=sigIn.yunits; 
                dataLabels=sigIn.dataLabels;
                plotProps=sigIn.plotProps;
                evalstring = strcat('sigOut=',class(sigIn),'(time, data,name, xlabelval, xunits, yunits,dataLabels,plotProps);');
                eval(evalstring);
                %sigOut = SignalObj(time, data,name, xlabelval, xunits, yunits,dataLabels,plotProps);
                sigOut.dataMask = sigIn.dataMask;
        end
        function sObjOut=resample(sObj, newSampleRate)
            %sObjOut=resample(sObj, newSampleRate)
            %sObjOut is a copy of sObj with the newSampleRate specified;
            if(sObj.sampleRate ~=newSampleRate)
                sObjOut = sObj.copySignal;
                if(or(~isnan(sObjOut.sampleRate),size(sObjOut.data,1)>1))
                    sObjOut.resampleMe(newSampleRate);
                end
            else
                sObjOut = sObj.copySignal;
            end
        end
        function resampleMe(sObj, newSampleRate)
            %resampleMe(sObj, newSampleRate)
            %Actually changes the sampleRate of the sObj signal. There is
            %no output variable. Recall that SignalObj's are handles.
            
            if(sObj.sampleRate~=newSampleRate)
                sObj.restoreToOriginal; %use default data for changing samplerate
                %Change the sampling rate of this object
                minTime=sObj.minTime;
                maxTime=sObj.maxTime;
                newTime=minTime:1/newSampleRate:maxTime;
                newData=zeros(length(newTime),sObj.dimension);
                for i=1:sObj.dimension
    %                 newData(:,i)= interp1(sObj.time,sObj.data(:,i),newTime,'spline','extrap');
                    newData(:,i)= interp1(sObj.time,sObj.data(:,i),newTime,'nearest',0);
                end
                sObj.time=newTime';
                sObj.data=newData;
                sObj.sampleRate=newSampleRate;
            end
        end
        function restoreToOriginal(sObj,rMask)
            %restoreToOriginal(sObj)
            %return the signal to its original state when created. The
            %dataMasks are the only thing preserved.
            %
            %restoreToOriginal(sObj,rMask)
            %if rMask=1, then the dataMask is reset. 
            
            if(nargin<2)
              rMask = 0; %keep mask even when data reset
            end
             [time,data]=sObj.getOriginalData;
             sObj.time=time;
             sObj.data=data;
             sObj.minTime=min(time);
             sObj.maxTime=max(time);
             sObj.sampleRate = 1/mean(diff(time));
             if(rMask==1)
                 sObj.resetMask;
             end
        end
        function resetMask(sObj)
            %resetMask(sObj)
            %Resets the dataMask for the SignalObj so that all components
            %of the object are visible.
            sObj.dataMask = ones(1,sObj.dimension);       
        end
        function ind = findIndFromDataMask(sObj)
            %ind = findIndFromDataMask(sObj) 
            %returns the indices of the visible components of the data
            ind = find(sObj.dataMask ==1);
%             if(isempty(ind))
%                 error('All data masked out of SignalObj');
%             end
        end        
        function ind = findNearestTimeIndices(sObj,times)
            ind = zeros(size(times));
            for i=1:length(ind)
                ind(i) = sObj.findNearestTimeIndex(times(i));
            end
        end
        
        function ind=findNearestTimeIndex(sObj, time)
            %ind=findNearestTimeIndex(sObj, time)
            %returns the nearest index to t=time;
            if(time<sObj.minTime)
                ind=1;
            elseif(time>sObj.maxTime)
                ind=length(sObj.time);
            else
                ind1=find(sObj.time>=time,1,'first');
                ind2=find(sObj.time<=time,1,'last');
                if(abs(sObj.time(ind1)-time)<=abs(sObj.time(ind2)-time))
                    ind=ind1;
                else
                    ind=ind2;
                end
%                 if(time-sObj.minTime < sObj.maxTime-time) %Time is close to beginning
%                     ind=min([ind1, ind2]);
%                 else
%                     ind=max([ind1, ind2]);
%                 end
            end
        end
        function sOut = shift(sObj, deltaT,updateLabels)
            %sOut = shift(sObj, deltaT
            %returns a SignalObj sOut shifted by deltaT. If deltaT is positive
            %the SignalObj is moved deltaT time steps forward
            %if deltaT is negative, the SignalObj is moved deltaT units
            %backwards.
            if(nargin<3)
                updateLabels =0;
            end
           
            %compute number of samples delayed
            sOut=sObj.copySignal;
            if(deltaT~=0)
                newMinTime = sOut.minTime+deltaT;
                newMaxTime = sOut.maxTime+deltaT;
                newTime=newMinTime:1/sOut.sampleRate:newMaxTime;
                sOut.time = newTime;
                sOut.minTime = newMinTime;
                sOut.maxTime = newMaxTime;
                if(updateLabels)
                    dataLabels = strcat(sObj.dataLabels,'(t-',num2str(deltaT),')');
                    sOut.setName(strcat(sObj.name,'(t-',num2str(deltaT),')'));
                    sOut.setDataLabels(dataLabels);
                end
                        
            end
        end
        function shiftMe(sObj,deltaT,updateLabels)
            %shiftMe(sObj,deltaT)
            % same as shift(sObj,deltaT) except that no signal is return.
            % The shift is done to sObj.
            if(nargin<3)
                updateLabels=0;
            end
            sTemp = sObj.shift(deltaT,updateLabels);
            sObj.data=sTemp.data;
            sObj.time=sTemp.time;
        end        
        function alignTime(sObj, timeMarker,newTime)
            if(sObj.minTime<=timeMarker && sObj.maxTime>=timeMarker)
                deltaT=newTime-timeMarker;
                sObj.shiftMe(deltaT);            
            end
        end     
        function answer = plotPropsSet(sObj)
            %answer = plotPropsSet(sObj)
            %returns 1 if any of the plotting properties has been set. Else
            %0.
            answer =0;
            for i=1:sObj.dimension
                if(~strcmp(sObj.getPlotProps(i),''))
                    answer=1;
                    break;
                end
            end
        end
        function answer = areDataLabelsEmpty(sObj)
            %answer = areDataLabelsEmpty(sObj)
            %returns 0 is dataLabels are empty and 1 otherwise
            answer = 1;
            for i=1:length(sObj.dataLabels);
                if(~strcmp(sObj.dataLabels{i},''))
                    answer = 0;
                    break;
                end
            end
        end
        function answer = isLabelPresent(sObj, label)
            %answer = isLabelPresent(sObj, label)
            %returns 1 if label is present in sObj, and 0 otherwise
            if(isa(label,'char'))
                if(strcmp(label,'all')||~isempty(sObj.getIndexFromLabel(label)))
                    answer=1;
                else
                    answer=0;
                end
            else
                error('Labels must be a char');
            end
        end
        function answer = isMaskSet(sObj)
            %answer = isMaskSet(sObj)
            %returns 0 is mask is not set, 1 if it is.
            answer=any(sObj.dataMask==0);
        end
        function sArray = convertNamesToIndices(sObj, selectorArray)
            %sArray = convertNamesToIndices(sObj, selectorArray)
            %converts the names in selectorArray to a vector of indices.
            
            if(sObj.areDataLabelsEmpty)
                sArray = 1:sObj.dimension; % Return all the data;
                fprintf('tried to find data by labels but data doesnot have labels assigned');
            else
                if(isa(selectorArray, 'char'))
                    if(strcmp(selectorArray,'all'))
                        sArray=1:sObj.dimension;
                    elseif(sObj.isLabelPresent(selectorArray))
                        sArray=sObj.getIndexFromLabel(selectorArray);
                    else
                        error('Specified label does not match data label');
                    end
                        
                elseif(isa(selectorArray, 'double'))
                    sArray=selectorArray;
                elseif(isa(selectorArray, 'cell'))
                    sArray = zeros(1, length(selectorArray));
                    for i=1:length(sArray)
                        if(sObj.isLabelPresent(selectorArray{i}))
                            sArray(i) = sObj.getIndexFromLabel(selectorArray{i});
                        end
                    end
                else
                    error('selectorArray cells must contain text');
                end 
            end
        end
        function [sAligned,meanTime] = alignToMax(sObj)
            
            [indices,values] = sObj.findGlobalPeak('maxima');
            meanTime = mean(indices);
            deltaT = -(indices-meanTime); % if index is greater than mean index then deltaT is negative so that we can align it witht the mean
            
            for i=1:sObj.dimension
               if(i==1)
                   sAligned = sObj.getSubSignal(i).shift(deltaT(i));
               else
                   sAligned = sAligned.merge(sObj.getSubSignal(i).shift(deltaT(i)));
               end
            end
        end
        function [ind, val] = findGlobalPeak(sObj,type)
           if(nargin<2)
               type='maxima';
           end            
           if(strcmp(type,'maxima'))
               [val,index] = max(sObj.data);
               ind = sObj.time(index);
           elseif(strcmp(type,'minima'))
               [val,index] = min(sOBj.data);
               ind = sObj.time(index);
           end
        end
        function [indices, values] = findPeaks(sObj,type,minDistance)
            %[indices, values] = findPeaks(sObj,type)
            %type:'minima' or 'maxima'
            %indices: indices at which the minima or maxima occur
            %values: values at the minima or maxima
            if(nargin<3)
                minDistance = round(sObj.sampleRate*(sObj.maxTime-sObj.minTime)/10);
            end
            if(nargin<2)
                type='maxima';
            end
            values=cell(1,sObj.dimension);
            indices=cell(1,sObj.dimension);
            if(strcmp(type,'maxima'))
                for i=1:sObj.dimension
                    [values{i},indices{i}] = findPeaks(sObj.data(:,i),'MINPEAKDISTANCE',minDistance);
                end
            elseif(strcmp(type,'minima'))
                for i=1:sObj.dimension
                    [values{i},indices{i}] = findPeaks(sObj.data(:,i),'MINPEAKDISTANCE',minDistance);
                end
            end
        end

        function [indices, values] = findMaxima(sObj)
            %Same as findPeaks(sObj,'maxima');
            [indices,values]=sObj.findPeaks('maxima');
        end
        function [indices, values] = findMinima(sObj)
            %Same as findPeaks(sObj,'minima');
            [indices,values]= sObj.findPeaks('minima');
        end        
        function clearPlotProps(sObj,index)
            %clearPlotProps(sObj,index)
            %clear the plotProps{index} if index is specified.
            %Otherwise all plotProps are cleared.
            if(nargin<2)
                index=1:sObj.dimension;
            end
            tempCell = cell(length(index),1);
            for i=index
                sObj.plotProps{i} = cell2str(tempCell{i});
            end
            
        end
%         function v=get.vars(sObj)
%             %returns a structure with fieldnames for each dataLabel
%             %eg: if s has components x and y. v will have two fields: v.x
%             %and v.y .
%             %Accessing v.x returns a SignalObj containing just the data
%             %corresponding to the label x.
%             v=[];
%             if(~sObj.areDataLabelsEmpty)
%                 UniqueSigLabels = unique(sObj.dataLabels);
%                 
%                 for i=1:length(UniqueSigLabels)
%                     if(~isempty(UniqueSigLabels{i}))
%                         K=strfind(UniqueSigLabels{i},'\');
%                         actString = UniqueSigLabels{i}(K+1:end);
%                         eval(strcat('v.(''',actString,''')=sObj.getSubSignal(''',UniqueSigLabels{i},''');'))
%                     end
%                         
%                 end
%             %sObj.vars=v;
%             else
%                 for i=1:sObj.dimension
%                     if(~isempty(sObj.dataLabels{i}))
%                         eval(strcat('v.(''',sObj.dataLabels{i},''')=sObj.getSubSignal(i);'))
%                     end
%                 end
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Change of Representation Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function structure=dataToStructure(sObj,selectorArray)
            %structure=dataToStructure(sObj)
            %return unmasked data in sObj as a standard matlab structure
            %with fields:
            %structure.time
            %structure.signals.values
            %
            %structure=dataToStructure(sObj,selectorArray)
            %returns a structure containing only the components specified
            %by the selectorArray
            
            if(nargin<2)
                if(sObj.isMaskSet)
                    selectorArray=sObj.findIndFromDataMask;
                else
                    selectorArray=1:sObj.dimension;
                end
            end
            %Convert to a standard matlab structure
            structure.time=sObj.time;
            structure.signals.values=sObj.dataToMatrix(selectorArray);
            structure.name        = sObj.name; 
            structure.dimension   = min(sObj.dimension,length(selectorArray));
            structure.signals.dimensions = min(sObj.dimension,length(selectorArray));
            structure.minTime     = sObj.minTime;
            structure.maxTime     = sObj.maxTime;
            structure.xlabelval   = sObj.xlabelval;
            structure.xunits      = sObj.xunits;
            structure.yunits      = sObj.yunits;
            structure.dataLabels  = sObj.dataLabels(selectorArray);
            structure.dataMask    = sObj.dataMask(selectorArray);
            structure.sampleRate  = sObj.sampleRate;
            structure.plotProps   = sObj.plotProps(selectorArray);
        end
        function dataMat=dataToMatrix(sObj,selectorArray)
            %dataMat=dataToMatrix(sObj)
            %returns the data for the visible components (i.e. those not
            %masked out, in a nxm matrix where n is the determined by the
            %sampleRate and the length of the current time vector, and m is
            %the number of unmasked dimensions.
            %
            %dataMat=dataToMatrix(sObj,selectorArray)
            %returns a matrix containing only the components specified
            %by the selectorArray
            if(nargin<2)
                if(sObj.isMaskSet)
                    selectorArray=sObj.findIndFromDataMask;
                else
                    selectorArray=1:sObj.dimension;
                end
            end
            % Convert data to a standard matrix. 
            % The matrix is length(times)x dimension 
            dataMat=sObj.data(:,selectorArray);
        end
        function sOut = getSubSignal(sObj,identifier)
            %sOut = getSubSignal(sObj,identifier)
            %Returns a new signal that consists of only the data
            %correspnding to the identifier. The identifier can be a string specifying the name
            %of the component, a cell array of strings specifying various
            %components, a vector specifying the indices of the components
            %of interest.
            if(isa(identifier,'cell'))
                if(isa(identifier{1},'char'))
                    sOut = sObj.getSubSignalFromNames(identifier);
                elseif(isa(identifier{1},'double'))
                    sOut = sObj.getSubSignalFromInd(identifier);
                else
                    error('Cells must contain strings!');
                end
            elseif(isa(identifier,'char'))
                sOut = sObj.getSubSignalFromNames(identifier);
            elseif(isa(identifier,'double'))
                sOut = sObj.getSubSignalFromInd(identifier);
            end
        end 
        function sOut = normWindowedSignal(sObj,windowTimes,numPoints,lbound,ubound)
            %Use the lower and upper bounds to specify the smallest and
            %largest window sizes that are allowed. Windows smaller or
            %larger than these will be ignored.
            if(nargin<5)
                if(nargin>4)
                    ubound = lbound;
                else
                    ubound = [];
                end
                
            end
            if(nargin<4 || isempty(lbound));
                lbound=[];
            end
            if(nargin<3 || isempty(numPoints))
                numPoints =100;
            end
            
          if(sObj.dimension == 1)
              data=[];
              for i=1:length(windowTimes)-1 
                  
                minTime = windowTimes(i);
                maxTime = windowTimes(i+1);
                if(and(~isempty(lbound),~isempty(ubound)))
                    if(and(abs(maxTime-minTime)<=(ubound),abs(maxTime-minTime)>=(lbound)))
                        dim=size(data,2);
                        time = linspace(minTime,maxTime,numPoints);
                        actSig= sObj.getSigInTimeWindow(minTime, maxTime);
%                         data(:,dim+1) = interp1(actSig.time,actSig.data,time,'spline')';

                        data(:,dim+1) = interp1(actSig.time,actSig.data,time,'nearest',0)';
                    end
                else
                   time = linspace(minTime,maxTime,numPoints);
                   actSig= sObj.getSigInTimeWindow(minTime, maxTime);
%                    data(:,i) = interp1(actSig.time,actSig.data,time,'spline')';
                    data(:,i) = interp1(actSig.time,actSig.data,time,'nearest',0)';
                end
              end  
                  actTime = (0:1:numPoints-1)./numPoints; %time -time(1);
                  name = sObj.name;
                  xlabelval = sObj.xlabelval;
                  xunits = '%';
                  yunits = sObj.yunits;
                  dataLabels(1:size(data,2)) =sObj.dataLabels;
                  %plotProps = sObj.plotProps;
%                   if(i==1)
                        evalstring = strcat('sOut=',class(sObj),'(actTime, data,name, xlabelval, xunits, yunits,dataLabels);');
                        eval(evalstring);
%                         sOut.setMinTime(0);
%                     else
%                         evalstring = strcat('sOut = sOut.merge(',class(sObj),'(actTime, data,name, xlabelval, xunits, yunits,dataLabels,plotProps));');
%                         eval(evalstring);
%                     end
                    
          end
        end
        
        function sOut = windowedSignal(sObj,windowTimes)
            if(sObj.dimension == 1)
                for i=1:length(windowTimes)-1
                    if(i==1)
                        sOut = sObj.getSigInTimeWindow(windowTimes(1), windowTimes(2));
                    else
                        temp =sObj.getSigInTimeWindow(windowTimes(i), windowTimes(i+1));
                        temp = temp.shift(-windowTimes(i));
                        sOut=sOut.merge(temp);
                    end
                end  
                
            end
        end
        function wSignals = getSigInTimeWindow(sObj,wMin,wMax,holdVals)
            %wSignal = getSigInTimeWindow(sObj)
            %returns a copy of sObj from t0=sObj.minTime to tf=sObj.maxTime;
            %wSignal = getSigInTimeWindow(sObj,wMin)
            %returns a new signal starting at t0=wMin to tf=sObj.maxTime
            %
            %wSignal = getSigInTimeWindow(sObj,wMin,wMax)
            %returns a new signal starting at t0=wMin to tf=wMax;
            %
            %wSignal = getSigInTimeWindow(sObj,wMin,wMax,holdVals)
            %if holdVals =1 and wMin or wMax are outside the bounds of the
            %current SignalObj time vector, then the endpoint values are
            %held. Otherwise, this region is padded with zeros;
            if(nargin<4)
                holdVals =0;
            end
            if(nargin<3)
                wMax=sObj.maxTime;
            end
            if(nargin<2)
                wMin=sObj.minTime;
            end
            
            if(length(wMin)~=length(wMax))
                error('Window minTimes must contain the same number of elements as window maxTimes');
            end
            if(length(wMin)==1 && sObj.minTime==wMin && sObj.maxTime==wMax)
                wSignals = sObj.copySignal;
            else
                for i=1:length(wMin)
                    wSignal = sObj.copySignal;
                    if(wMin(i)<wSignal.minTime) 
                        %if the window starts before the SignalObj pad with zeros
                        %and set as new starttime
                        wSignal.setMinTime(wMin(i),holdVals);
                    end
                    if(wMax(i)>wSignal.maxTime)
                        %if wMax is too big, pad with zeros at end and set as new
                        %end
                        wSignal.setMaxTime(wMax(i),holdVals); 
                    end

                    startIndex= wSignal.findNearestTimeIndex(wMin(i));
                    endIndex  = wSignal.findNearestTimeIndex(wMax(i));
                    wSignal.time=wSignal.time(startIndex:endIndex); 
                    wSignal.data=wSignal.data(startIndex:endIndex,:);
                    if(length(wMin)>1)
                        for j=1:wSignal.dimension
                            dataLabels{j} = strcat(wSignal.dataLabels{j},'_{',num2str(i),'}');
                        end
                    else
                        if(wSignal.dimension==1)
                            dataLabels = wSignal.dataLabels;
                        else
                            for j=1:wSignal.dimension
                                dataLabels{j} = wSignal.dataLabels{j};
                            end
                        end
                    end
                    wSignal.setDataLabels(dataLabels);
                    wSignal.setMinTime;
                    wSignal.setMaxTime;
                    if(i==1)
                        wSignals = wSignal;
                    else
                        wSignals = wSignals.merge(wSignal);
                    end
                end 
            end
        end        

        function [s sigIndex]=getSubSignalsWithinNStd(sObj,nStd)
            mSig=mean(sObj.data,2);
            stdSig = std(sObj.data,0,2);
            minVal = mSig-nStd*stdSig;
            maxVal = mSig+nStd*stdSig;
            indMin= find(sum(sObj.data-minVal*ones(1,sObj.dimension)<0)==0);
            indMax=find(sum(sObj.data-maxVal*ones(1,sObj.dimension)>0)==0);
            sigIndex = intersect(indMin,indMax);
            s=sObj.getSubSignal(sigIndex);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function h=plot(sObj,selectorArray,plotPropsIn,handle)
            
        %   h=plot(sObj,selectorArray,plotProps,handle)
        %   plotProps: is a cell array with strings in each a entry. 
        %              The strings must evaluate to parameters that can be 
        %              passed to the standard matlab plot function.
        %   example: To plot a three dimensional signals with each line
        %   with specific properties (not the use of double quotes on the
        %   color entries so that the actual arguments have a quotes around
        %   them. 
        %
        % >> plotProps = {{' ''-.g'', ''LineWidth'' ,2'},...
        %                 {' ''k'', ''LineWidth'' ,2'},...
        %                 {' ''-.b'' '}};
        %
        %   selectorArrray: is either an numeric array indexing the
        %                   elements within the SignalObj data or a cell array
        %                   of characters that refers to the data elements
        %                   by their labels.
        %   handle: is figure handle to tell the signal where to plot
        %           itself.
        % % 8-25-09 : Legend wont plot is SignalObj dimension >10;
     
           
           
           if(nargin<4)
              handle=gca;
           end
           if((nargin<3) || isempty(plotPropsIn))
              %plotPropsIn=[];
              plotPropsIn=sObj.plotProps;
           end
           if((nargin<2) || isempty(selectorArray))
                if(sObj.isMaskSet)
                    selectorArray=sObj.findIndFromDataMask;
                    if(isempty(selectorArray))
                        h=[];
                        return;
                    end                        
                else
                    selectorArray=1:sObj.dimension;
                end
           end


           
            if(isa(selectorArray,'cell'))
                sArray = sObj.convertNamesToIndices(selectorArray);
            elseif(isa(selectorArray,'char'))
                sArray = sObj.convertNamesToIndices(selectorArray);
            elseif(isa(selectorArray,'double'))
                sArray = selectorArray;
            else
                error('selectorArray must contain either dataLabels or their indices');
            end
            
            
           % if plot params are passed in then assign them to the SignalObj
           if(length(plotPropsIn)==length(sArray))
                for i=1:length(sArray)
                    sObj.setPlotProps(plotPropsIn{i},sArray(i));
                end
           elseif(length(plotPropsIn)<=length(sArray) && length(plotPropsIn)==1)
               for i=1:length(sArray)
                    sObj.setPlotProps(plotPropsIn,sArray(i));
               end
           end
           
           
              
          
            %If the dimensions of the data are larger than one can specify
            %the plot style for each dimension by placing their values in
            %an cell array.
%             h=gcf;figure(h);
            if((nargin>=3) && ~isempty(plotPropsIn))  %is we got called with parameter use those
                %sObj.setPlotProps(plotProps); %Only accepts entire cell of plotting properties as input
                %if(length(plotCell)>=length(sArray)) %if we specified plot parameters for each dimension
                for i=sArray
                    %plotStr=cell2str(plotCell{:,find(sArray==i)});
                    plotStr=cell2str(sObj.getPlotProps(i));
                    if(~strcmp(plotStr,'') && ~isempty(plotStr))
                        evalstring = strcat('h(',num2str(find(sArray==i)),')=plot(handle,sObj.time,sObj.data(:,',num2str(i),'),', plotStr,');');
                        eval(evalstring); hold on;
                    else
%                         set(gca,handle);
                        axes(handle);
                        evalstring = strcat('h(',num2str(find(sArray==i)),')=plot(sObj.time,sObj.data(:,',num2str(i),'));');
%                         evalstring = strcat('h(',num2str(find(sArray==i)),')=plot(handle,sObj.time,sObj.data(:,',num2str(i),'));');
                        eval(evalstring); hold all;
                    end
                    %evalstring
                    %eval(evalstring); hold on;
                end
            elseif(sObj.plotPropsSet) %% We use the values that have been set previously
                %plotCell = cell2str(sObj.getPlotProps(i));
                %fprintf('varargin is not empty')
                for i=sArray
                    plotStr =cell2str(sObj.getPlotProps(i));
                    %plotStr=cell2str(plotCell{:,1});
                    if(~strcmp(plotStr,'') && ~isempty(plotStr))
                        evalstring = strcat('h(',num2str(i),')=plot(handle,sObj.time,sObj.data(:,',num2str(i),'),', plotStr,');');
                    else
                        evalstring = strcat('h(',num2str(i),')=plot(handle,sObj.time,sObj.data(:,',num2str(i),'));');
                    end
                    %evalstring
                    if(~isempty(sObj.data(:,i)))
                        eval(evalstring); hold on;
                    end
                end
            else %We didnt get any plotting properties. Use matlab default
                %plotCell{1} = '''';
                %fprintf('varargin is empty')
                set(gcf,'CurrentAxes',handle);
                h=plot(handle,sObj.time,sObj.data(:,sArray));
            end

                sObj.setupPlots(handle,sArray);
        end        
        function setupPlots(sObj,handle,sArray)
            %setupPlots(sObj,sArray)
            %Sets the labels for the x-axis, y-axis, and legend.
            %sArray is an array on indexes corresponding to which labels
            %will appear in the legend.
%              set(gcf,'CurrentAxes',handle);
            warning off;
            if(~strcmp(sObj.xunits,''))
                xunitsStr=strcat('\; [',sObj.xunits,']'); %\; is a large space in latex
            else
                xunitsStr='';
            end
            
            if(~strcmp(sObj.yunits,''))
                yunitsStr=strcat('\; [',sObj.yunits,']');
            else
                yunitsStr='';
            end
            
            if(strcmp(sObj.xlabelval,''))
                xlabel(strcat('$$',xunitsStr,'$$'),'Interpreter','latex');
                
            else                
%                 strcat(['$$' sObj.xlabelval xunitsStr '$$'])
                xlabel(strcat('$$',[sObj.xlabelval xunitsStr],'$$'),'Interpreter','latex');
                
            end
            if(~strcmp(sObj.name,''))
                if(~strcmp(yunitsStr,''))
                    %strcat('$$',sObj.name,yunitsStr,'$$')
                    ylabel(strcat('$$',[sObj.name yunitsStr],'$$'),'Interpreter','latex');
                    
                else
                     
                    ylabel(strcat('$$',sObj.name,'$$'),'Interpreter','latex');
                    %ylabel(sObj.name,'Interpreter','none');
                    
                end
            end
            
            if(~sObj.areDataLabelsEmpty) %%&& sObj.dimension<10)
%                 labelArray= cell(1,length(sArray));
                labelArray =sObj.dataLabels(sArray);
%                 for i=1:length(sArray)
%                     labelArray{i} = strcat('$$',sObj.dataLabels{sArray(i)},'$$');
%                 end
                legend(handle,labelArray);%,'Interpreter','latex');
            end
            axis tight;
            warning on;
        end
        
        function h=plotVariability(sObj,selectorArray)
            %h=plotVariability(sObj)
            %Calls plotAllVariability(sObj) for each set of unique labels
            %in the SignalObj. For example, if a SignalObj s has for
            %components (2 with label 'x', and 2 with label 'y'), two
            %figures will be generated, one for the variability of the x
            %component and one for the variability of the y component. If
            %all components have the same label, then only one plot will be
            %generated. Plot handles are returned in h.
            %
            %h=plotVariability(sObj,selectorArray)
            %if selectorArray is specified, then plotAllVariability is
            %called once on a SignalObj that only has the components
            %specified by the selectorArray.
            if(nargin<2)
                if(~sObj.areDataLabelsEmpty)
                    uLabels =unique(sObj.dataLabels);
                    for i=1:length(uLabels)
                        selectorArray{i} = sObj.getIndicesFromLabels(uLabels{i});
                    end
                else
                    selectorArray = 1:sObj.dimension;
                end
            end
            if(isa(selectorArray,'cell')) %More than component has the same 
                [numSubSignals]=length(selectorArray);
                for i=1:numSubSignals
                    %figure;
                    h(i)=sObj.getSubSignal(selectorArray{i}).plotAllVariability(getAvailableColor(i));
                end
            elseif(isa(selectorArray,'double'))
                h=sObj.getSubSignal(selectorArray).plotAllVariability;
            end
        end
        
        function h=plotAllVariability(sObj,faceColor,linewidth,ciUpper,ciLower)
            %h=plotAllVariability(sObj)
            % Computes the mean across all components in the signal and
            % plots +/- 1 standard deviation from the mean. All other
            % parameters are optional.
            %
            % faceColor: defaults to red if not specified. 
            % linewidth: defaults to 3 if not specified. 
            %
            % ciUpper: if a single number then specifies the multiple of
            %          standard deviations to be used. i.e. ciUpper=1
            %          results in default. If ciUpper is an double array of
            %          the same length as the signal, then the upper bound
            %          of the confidence interval at time t will be
            %          ciUpper+mean(sObj). If ciUpper is a SignalObj, then
            %          the behavior is the same as if a double array were
            %          specified. If only ciUpper is specified, ciLower is
            %          assigned the same maginitude
            % ciLower: if a single number then specifies the multiple of
            %          standard deviations to be used. i.e. ciLower=1
            %          results in default. If ciLower is an double array of
            %          the same length as the signal, then the lower bound
            %          of the confidence interval at time t will be
            %          mean(sObj)-ciLower. If ciLower is a SignalObj, then
            %          the behavior is the same as if a double array were
            %          specified.
            if(nargin<4 || isempty(ciUpper))
                ciUpper =1.96;
            end
            if(nargin<=4)                
                ciLower=ciUpper;
            end   
            if(nargin<2 || isempty(faceColor))
                faceColor=getAvailableColor(1);
            end
            if(nargin<3 || isempty(linewidth))
                linewidth=3;
            end
%              sObj=s;
            meanSig = mean(sObj,2);
            stdSig  = std(sObj,[],2);
            
            if(length(ciUpper)==1)
                ciTop=meanSig+ciUpper.*stdSig;
            elseif(length(ciUpper)==length(sObj.time))
                ciTop=meanSig+ciUpper;
            else
                error('Upper confidence interval must be either length 1 or same length as the time vector');
            end
            
            if(length(ciLower)==1)
                ciBottom=meanSig-ciLower.*stdSig;
            elseif(length(ciUpper)==length(sObj.time))
                ciBottom=meanSig-ciLower;
            else
                error('Lower confidence interval must be either length 1 or same length as the time vector');
            end 
            
            ci2=ciTop.dataToMatrix;           
            ci1=ciBottom.dataToMatrix;
            

           
            p=patch([reshape(sObj.time,[1 length(sObj.time)]) reshape(fliplr(sObj.time'),[1 length(sObj.time)])],...
                    [reshape(ci1',[1 length(sObj.time)])       reshape(fliplr(ci2'),[1,length(sObj.time)])],strcat('',faceColor,''));
            set(p,'facecolor',faceColor,'edgecolor','none');
            alpha(.5);hold on;
            h=plot(sObj.time,meanSig.dataToMatrix,'k-','linewidth',linewidth);
            
            sArray=1; %only want the first label to plot;
            meanSig.setupPlots(sArray);
%             legend([h,p],'mean','variance');
            
        end
        

        
        
    end % public methods
    methods (Access = private)
        function index      = getIndexFromLabel(sObj,label)
            %index      = getIndexFromLabel(sObj,label)
            %Given a label, returns the index to the SignalObj component
            %with the specified label
            index=[];
            for i=1:length(sObj.dataLabels)
                if(strcmp(label, sObj.dataLabels{i}))
                    index = [index i];
                end
            end
        end 
    
        
        function setDataMask(sObj, dataMask)
        % setDataMask(sObj, dataMask)
        % dataMask is a vector of ones and zeros of length sObj.dimension
        % components marked by ones will be visible.i
        %DataMasks affect how the SignalObj is visualized and converted to
        %other representations. It doesnt change the dimensions of the
        %object, nor does it delete data. 
            if(length(dataMask)==sObj.dimension)
                  sObj.dataMask = dataMask;
            end
                
        end
        function setMaskByInd(sObj,index)
            % setMaskByInd(sObj,index)
            % if length(index)== sObj.dimension then the mask is set using
            % the values in index.
            % if length(index)< sObj.dimension only the component specified
            % by the index will be visible. All others will be hidden.
            if(length(index)==sObj.dimension)
                sObj.setDataMask(index);
            else
                mask=zeros(1,sObj.dimension);
                mask(index)=1;
                sObj.setDataMask(mask);
            end
        end
        function setMaskByLabels(sObj,labels)
            %setMaskByLabels(sObj,labels)
            %labels is a cell array with the labels of the components that
            %are to remain visible. The indices for these labels are found
            %and the mask set using setDataMask
            ind=sObj.getIndicesFromLabels(labels);
            mask=zeros(1,sObj.dimension);
            mask(ind)=1;
            sObj.setDataMask(mask);
        end      
    
        function sOut = getSubSignalFromInd(sObj, selectorArray)
            if(nargin<2)
                if(sObj.isDataMaskSet)
                    selectorArray=sObj.findIndFromDataMask;
                else
                    selectorArray=1:sObj.dimension;
                end
            end
            s3   = sObj.copySignal;
            time = s3.time;
            if(isa(selectorArray,'cell'))
                data=[];
               for i=1:length(selectorArray)
                   offset=0;
                   data = [data s3.data(:,selectorArray{i})];
                   for j=1:length(selectorArray{i})
                        dataLabels{offset + j} = s3.dataLabels{selectorArray{i}(j)};
                        if(~isempty(s3.plotProps))
                            plotProps{offset+ j} = s3.plotProps{selectorArray{i}(j)};
                        end
                   end
                   offset = offset + length(selectorArray{i});
               end
            else
                data = s3.data(:,selectorArray);
                dataLabels = cell(1,length(selectorArray));
                plotProps  = [];
                for i=1:length(selectorArray)
                    dataLabels{i}= s3.dataLabels{selectorArray(i)};
                    if(~isempty(s3.plotProps))
                        plotProps{i} = s3.plotProps{selectorArray(i)};
                    end
                end
            end
            
            
            

            name = s3.name;
            xlabelval = s3.xlabelval; xunits=s3.xunits; yunits = s3.yunits;
            
            %sOut=SignalObj(time, data, name, xlabelval, xunits, yunits, dataLabels, plotProps);
            evalstring = strcat('sOut=',class(sObj),'(time, data,name, xlabelval, xunits, yunits,dataLabels,plotProps);');
            eval(evalstring);
        end
        function sOut = getSubSignalFromNames(sObj,labels)
            ind   = sObj.getIndicesFromLabels(labels);
            if(isa(ind,'cell'))
                if(length(ind)>1)
                    for i=1:length(ind)
                        s{i} = getSubSignalFromInd(ind{i});
                    end
                    sOut = s{1}.merge(s{2:end});
                else
                    sOut = sObj.getSubSignalFromInd(ind{1});
                end
            elseif(isa(ind,'double'))
                sOut  = sObj.getSubSignalFromInd(ind);
            end
        end
    end %private methods
    methods (Static)
        function sObj =signalFromStruct(structure)
            fNames = {'time','signals','name','xlabelval','xunits','yunits','dataLabels','plotProps'};
            structField = fields(structure);
            for i=1:length(fNames)
                if(~any(strcmp(structField,fNames{i})))
                    error(['Field ' fNames{i} 'not present is structure! Needed to make a signal from a structure']);
                end
            end
            sObj=SignalObj(structure.time, structure.signals.values, structure.name, structure.xlabelval, structure.xunits, structure.yunits, structure.dataLabels, structure.plotProps);
            sObj.setMask(structure.dataMask);
        end 
        function simpleStructure = convertSigStructureToStructure(sigStructure)
           if(isa(sigStructure,'struct'))
                fNames = fields(sigStructure);
                for i=1:length(fNames)
                    if(isa(sigStructure.(fNames{i}),'SignalObj'))
                        simpleStructure.(fNames{i}) = sigStructure.(fNames{i}).dataToStructure;
                    end
                end
           end
        end
        function sigStructure = convertSimpleStructureToSigStructure(simpleStructure)
            if(isa(simpleStructure,'struct'))
                fNames = fields(simpleStructure);
                for i=1:length(fNames)
                    if(isa(simpleStructure.(fNames{i}),'struct'))
                        sigStructure.(fNames{i}) = SignalObj.signalFromStruct(simpleStructure.(fNames{i}));
                    end
                end
           end
        end
            
        
        
    end
end%classdef

% Helper function
function stringout=cell2str(input)
    if(isempty(input))
        stringout = '';
    else
        if(isa(input,'char'))
            stringout=input;
        elseif(isa(input,'cell'))
            stringout=cell2str(input{1});
        end
    end
end
function color=getAvailableColor(index)
    availableColors = {'g','b','r','y','c','m','k'};
    color=availableColors{mod(index,length(availableColors))+1};
end


