classdef History <handle
%HISTORY defines windows of interest for analysis. Given a set of
%windowTimes of length N, N-1 windows are created, these windows are:
% $w_i$ is a window from windowTimes(i) to windowTimes(i+1)
%
% Usage:
% HistObj = History(windowTimes)
%           Window times is a vector of times; We make sure that we have
%           them in order and use them to specify windows in which we are
%           interested in computing the history;
%        
% <a href="matlab: methods('History')">methods</a>
% <a href="matlab:web('HistoryExamples.html', '-helpbrowser')">History Examples</a> 
% Reference page in Help browser
% <a href="matlab: doc('History')">doc History</a>


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
        windowTimes % times specifying the desired history windows
        minTime 
        maxTime
    end
    
    methods 
          function HistObj = History(windowTimes,minTime,maxTime)
            % HistObj = History(windowTimes)
            % Window times is a vector of times; We make sure that we have
            % them in order and use them to specify windows in which we are
            % interested in computing the history;
            if(nargin<3)
                maxTime=[];
            end
            if(nargin<2)
                minTime=[];
            end
                
            HistObj.windowTimes = sort(windowTimes);
            HistObj.minTime     = minTime;
            HistObj.maxTime     = maxTime;
          end

          function cov = computeHistory(HistObj, nst,historyNum,tn)
              % cov = computeHistory(HistObj, nst)
              %
              % returns a CovColl if more than one neural spike train is
              % received as input. 
              %
              % nst can be a nstColl, a cell array of nspikeTrains or a
              % single nspikeTrain object.
              %
              % each output covariate will have the same number of dimensions as
              % the number of history windows. The firing history
              % corresponding to the window $w_i$ is the ith component of
              % the covariate
              % If tn is specified only compute the history at time tn
              
              if(nargin<4)
                  tn=[];
              end
              if(nargin<3)
                  historyNum=[];
              end
              
              if(isa(nst,'nstColl'))
                  temp = cell(1,nst.numSpikeTrains);
                  for i=1:nst.numSpikeTrains
%                       if(strcmp(nst.getNST(i).name,''))
%                           nst.getNST(i).setName(strcat('n',num2str(i)));
%                       end
                      temp{i} =HistObj.computeNSTHistoryWindow(nst.getNST(i),historyNum,tn);
                      if(strcmp(temp{i}.name,'History'))
                          %then nspikeTrain didnt have a name number
                          %history for the covariate collection
                          temp{i}.setName(['History #' num2str(historyNum) ' for ' nst{i}.name]);
                      end
                  end
                  cov=CovColl(temp);
              elseif(isa(nst,'cell') && isa(nst{1},'nspikeTrain')) % a cell collection of neural spike trains
                  temp = cell(1,length(nst));
                  for i=1:length(nst)
                      temp{i} = HistObj.computeNSTHistoryWindow(nst{i},historyNum,tn);
                      if(strcmp(temp{i}.name,'History'))
                          %then nspikeTrain didnt have a name number
                          %history for the covariate collection
                          temp{i}.setName(['History #' num2str(historyNum) ' for ' nst{i}.name]);
                      end
                  end
                  cov = CovColl(temp);
                  cov.setSampleRate(nst{1}.sampleRate);
              elseif(isa(nst,'nspikeTrain'))
                  temp=HistObj.computeNSTHistoryWindow(nst,historyNum,tn);
                  temp.setName(['History #' num2str(historyNum) ' for ' nst.name]);
                  cov = CovColl(temp);
              else
                  error('Can only compute History for nstColl, cells, or nspikeTrain');
              end
              %size(cov.standardRep.data)
              if(~isa(nst,'cell'))
                cov.setSampleRate(nst.sampleRate);  
              end
          end
 
        function setWindow(HistObj,windowTimes)
            % setWindow(HistObj,windowTimes)
            % replaces HistObj.windowTimes with the windowTimes vector that
            % is being specified.
            HistObj.windowTimes = windowTimes;
        end
        
        function plot(HistObj)
            % plots each of the history windows
            tmin=HistObj.windowTimes(1:end-1);
            tmax=HistObj.windowTimes(2:end);
            sampleRate = 1000;
            data=zeros((max(tmax)-min(tmin))*sampleRate,length(tmax));
            for i = 1:length(tmax)
                indMin = max(1,(tmin(i)-min(tmin))*sampleRate);
                indMax = (tmax(i)-min(tmin))*sampleRate;
                data(indMin:indMax,i)=1;
                dataLabels{i} = strcat('[',num2str(tmin(i),3),',',num2str(tmax(i),3),']');
            end
                name='History';
                time=linspace(min(tmin),max(tmax),length(data));
                xlabelval = 'time'; xunits='s'; yunits='';
                s = SignalObj(time,data,name,xlabelval, xunits, yunits, dataLabels);
                s.plot; hold on;
                
                
        end
        function filterMat = toFilter(HistObj,delta)
%             if(nargin<2)
%                 delta = .001;
%             end
            

            tmin=HistObj.windowTimes(1:end-1);
            tmax=HistObj.windowTimes(2:end);
            timeVec=min(tmin):delta:max(tmax);
            
            a=ones(length(timeVec),1);
            b=zeros(length(tmax),length(timeVec));
            
            for i=1:length(tmax)
            
                
                NumSamples = ceil(tmax(i)/delta);
%               
                StartSample = ceil(tmin(i)/delta) +1; 
                

                b(i,(StartSample:NumSamples)+1)=1; %delay by 1
                den{i,1} = a(i);
                num{i,1} = b(i,:);
                
            end
            filterMat = tf(num,den,delta,'Variable','z^-1');
            
        end
        function structure = toStructure(HistObj)
            fNames = fieldnames(HistObj);
            for i=1:length(fNames)
               structure.(fNames{i}) = HistObj.(fNames{i}); 
            end
            
        end
        
    end
    methods (Static)
        function HistObj = fromStructure(structure)
            if(~isempty(structure))
                fNames = fieldnames(structure);
                windowTimes = structure.windowTimes;
                minTime = structure.minTime;
                maxTime = structure.maxTime;
                HistObj = History(windowTimes,minTime,maxTime);
            else
               HistObj = []; 
            end
        end
    end
    methods (Access = private)
         function cov = computeNSTHistoryWindow(HistObj,nst,historyNum,tn)
             if(nargin<4)
                 tn=[];
             end
             if(nargin<3)
                 historyNum=[];
             end
             
              s = nst.getSigRep;
              tmin=HistObj.windowTimes(1:end-1);
              tmax=HistObj.windowTimes(2:end);
              %get signal representionat from nst
              %find number of samples to first point in history window
              %find number of sample to second....
              data=[];
              dataLabels=cell(1,length(tmax));
              for i=1:length(tmax)
                  a=1;
                  
%                 b=zeros(1,s.findNearestTimeIndex(tmax(i)));
                  NumSamples = ceil(tmax(i)*nst.sampleRate);
                  b=zeros(1,NumSamples);
                  StartSample = ceil(tmin(i)*nst.sampleRate) +1; 
%                   b(s.findNearestTimeIndex(tmin(i)):(s.findNearestTimeIndex(tmax(i))-1))=1;

                  b(StartSample:NumSamples)=1;
                  sTemp = s.filter(b,a);
                  %Delay by 1 to make lag the actual spike 
                  bdelay=[0 1]; adelay=1;
                  sOut{i} = sTemp.filter(bdelay,adelay);
                  data=[data, sOut{i}.dataToMatrix];
                  if(isempty(historyNum))
                    dataLabels{i} = strcat('[',num2str(tmin(i),3),',',num2str(tmax(i),3),']');
                  else
                    dataLabels{i} = strcat('[',num2str(tmin(i),3),',',num2str(tmax(i),3),']_',num2str(historyNum));
                  end
              end
              %name =['History \; for \;', s.name]; %w\; is  a thick space in latex
              if(isempty(nst.name))
                name = 'History';
              else
                  name =['History ' nst.name];
              end
              xlabelval = s.xlabelval;
              xunits = s.xunits;
              yunits = s.yunits;
        
         
              if(~isempty(data))
                  if(~isempty(tn))
                      
                      cov = Covariate(s.time, data, name, xlabelval,xunits,yunits,dataLabels);
                      dataVals = cov.getValueAt(tn);
%                       figure(2); cov.plot; hold on; plot(tn-1/s.sampleRate, dataVals,'o'); pause(.001);
%                       cov = Covariate(tn-1/s.sampleRate, dataVals, name, xlabelval,xunits,yunits,dataLabels);
%                       
                      
                  else
                      cov = Covariate(s.time, data, name, xlabelval,xunits,yunits,dataLabels);
                      %Window the data if minTime and maxTime have been set.
                      minTime=HistObj.minTime;
                      maxTime=HistObj.maxTime;
                      if(isempty(minTime))
                          minTime=cov.minTime;
                      end
                      if(isempty(maxTime))
                          maxTime=cov.maxTime;
                      end
                      cov.resample(nst.sampleRate);
                      wCov=cov.getSigInTimeWindow(minTime,maxTime);
                      wCov.setMinTime(nst.minTime);
                      wCov.setMaxTime(nst.maxTime);
                      cov=wCov;
                  end
            
                  
              else
                  cov = Covariate([], data, name, xlabelval,xunits,yunits,dataLabels);
              end
              
             
              
             
              
              
              %size(data)
          end
    end
    
    
end

