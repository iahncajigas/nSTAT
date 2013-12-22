classdef CIF < handle
    %CIF - Conditional Intensity function. 
    %<a href="matlab: methods('CIF')">methods</a>
    % 
    %Reference page in Help browser
    %<a href="matlab: doc('CIF')">doc CIF</a>

    
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
        
        b           %Regression Coefficients
        varIn       %The labels for the coefficients in b
        stimVars    %The subset of labels that correspond to the stimulus
        indepVars   %
        stats
        fitType     % binomial or poisson - determines how lambda is related to the regression coefficients
        lambdaDelta % symbolic expression for the product of lambda and delta
        lambdaDeltaGamma
        LogLambdaDeltaGamma
        spikeTrain
        
        gradientLambdaDelta
        gradientLogLambdaDelta %symbolic expression for first partial w.r.t. to stimulus variables
        gradientLambdaDeltaGamma %symbolic expression for first partial w.r.t. to history coefficient variables
        gradientLogLambdaDeltaGamma %symbolic expression for first partial w.r.t. to history coefficient variables
        
        jacobianLambdaDelta
        jacobianLogLambdaDelta %symbolic expression for second partial w.r.t. to stimulus variables
        jacobianLambdaDeltaGamma %symbolic expression for second partial w.r.t. to history variables        
        jacobianLogLambdaDeltaGamma %symbolic expression for second partial w.r.t. to history variables
        history
        histCoeffs
        histCoeffVars %Defined in case we want to take derivatives with respect to the history params (as in the M-step of EM)
        histVars
        historyMat
%     end
    
%     properties (Hidden)

        lambdaDeltaFunction % function handle to evaluate lambda*delta
        lambdaDeltaGammaFunction
        LogLambdaDeltaGammaFunction
        gradientFunction % partial derivative of log(lambda*delta) w.r.t stimulus variables
        gradientLogFunction
        gradientFunctionGamma % partial derivative of lambda*delta w.r.t stimulus variables
        gradientLogFunctionGamma % partial derivative of log(lambda*delta) w.r.t stimulus variables
        
        jacobianFunction % second partial derivative of (lambda*delta) w.r.t. to stimulus variables
        jacobianLogFunction
        jacobianFunctionGamma % second partial derivative of lambda*delta w.r.t. to stimulus variables
        jacobianLogFunctionGamma % second partial derivative of log(lambda*delta) w.r.t. to stimulus variables
        
        argstr % parse out stimulus variables by each element since the above functions dont take vector inputs
        argstrLDGamma
        
    end
    
    methods
        function cifObj = CIF(beta,Xnames,stimNames,fitType,histCoeffs,historyObj,nst)
            % cifObj = CIF(beta,Xnames,stimNames,fitType)
            % beta: regression coefficients
            % 
            % Xnames: names of the variables in the order they are
            %         specified by beta.
            % 
            % stimNames: names of the subset of variables that are define
            %            the stimulus.
            % 
            % fitType: poisson or binomial - defines how the parameters are
            %          related to the CIF. For poisson, lamda*delta =
            %          exp(X*beta). For binomial, lambda*delta=
            %          exp(X*beta)/(1+exp(X*beta));
            %
            % histCoeffs: coefficients for each of the history windows
            %             defined in historyObj
            %
            % historyObj: an object of class History that defines the how
            %             the spiking activity is being windowed. This
            %             input could also be a vector of windowTimes to be
            %             used in creating the historyObj.
            
            if(nargin<7)
               cifObj.spikeTrain = []; 
            else
               cifObj.spikeTrain = nst.nstCopy;
            end
            if(nargin<6)
                cifObj.history=[];
            else
                cifObj.setHistory(historyObj);
            end
            if(nargin<5)
                cifObj.histCoeffs = [];
            else
                [r,c] = size(histCoeffs);
                if(r==1)
                    cifObj.histCoeffs = histCoeffs;
                elseif(c==1)
                    cifObj.histCoeffs = histCoeffs';
                else 
                    error('History Coefficient vector must have one dimension equal to 1');
                end
            end
            
            if(nargin<4)
                fitType = 'poisson';
            end
            
            
            if(isa(Xnames,'sym'))
                XnamesTemp=cell(length(Xnames),1);
                for i=1:length(beta)
                    XnamesTemp{i} = char(Xnames(i));
                end    
                Xnames=XnamesTemp;
            end
            
            % Define input variables as a vector;
            [r,c] = size(Xnames);
            if(r==1)
                Xnames = Xnames';
                cifObj.varIn = sym(Xnames);
            elseif(c==1)
                cifObj.varIn = sym(Xnames);
                
            else
                error(' Must have one dimension equal to 1');
                
            end
            
            % Define stimulus variables as a vector
            [r,c] = size(stimNames);
            if(r==1)
                cifObj.stimVars = sym(stimNames');
            elseif(c==1)
                cifObj.stimVars = sym(stimNames);
                
            else
                error(' Must have one dimension equal to 1');
                
            end
            
            % Define beta as a row vector
            if(isnumeric(beta))
                [r,c] = size(beta);
                if(r==1)
                    cifObj.b = beta;
                elseif(c==1)
                    cifObj.b = beta';
                elseif(isempty(beta))
%                     error('Coefficient vector must have one dimension equal to 1');
                      %define beta as symbolic
                      betaLabel=cell(1,length(cifObj.varIn));
                      for i=1:length(cifObj.varIn)
                          betaLabel{i} = strcat('b',num2str(i));
                      end
                      display('Beta is being treated as symbolic! Must provide an input vector length(beta)+length(Xnames) to evaluate');
                      
                      cifObj.b = sym(betaLabel);
                      beta = cifObj.b;
                      allVarNames = cell(length(Xnames)+length(betaLabel),1);
                      allVarNames(1:length(betaLabel)) = betaLabel;
                      allVarNames((length(betaLabel)+1):(length(betaLabel)+length(Xnames)))=Xnames;
                      cifObj.varIn =  sym(allVarNames);
                end
            elseif(isa(beta,'cell'))
                [r,c] = size(beta);
                if(r==1)
                    betaLabel = beta;
                    
                elseif(c==1)
                    betaLabel = beta';
                
                else
                    error(' Beta Must have one dimension equal to 1');
                end
                    cifObj.b = sym(betaLabel);
                    beta = cifObj.b;
                    allVarNames = cell(length(Xnames)+length(betaLabel),1);
                    allVarNames(1:length(betaLabel)) = betaLabel;
                    allVarNames((length(betaLabel)+1):(length(betaLabel)+length(Xnames)))=Xnames;
                    cifObj.varIn =  sym(allVarNames);
            elseif(isa(beta,'sym'))
                betaLabel=cell(1,length(beta));
                for i=1:length(beta)
                    betaLabel{i} = char(beta(i));
                end
                cifObj.b = sym(betaLabel);
                beta = cifObj.b;
                allVarNames = cell(length(Xnames)+length(betaLabel),1);
                allVarNames(1:length(betaLabel)) = betaLabel;
                allVarNames((length(betaLabel)+1):(length(betaLabel)+length(Xnames)))=Xnames;
                cifObj.varIn =  sym(allVarNames);
            end
                
            
            %Define History variables if they were passed in
            if(and(~isempty(cifObj.histCoeffs),~isempty(cifObj.history)))
                for i=1:length(cifObj.histCoeffs)
                    histNames{i} = strcat('dN',num2str(i));
                    histCoeffVars{i} = strcat('gamma',num2str(i));
                end
                cifObj.histVars = sym(histNames');
                cifObj.histCoeffVars = sym(histCoeffVars);
                histCoeffsVarsTrans = sym(histCoeffVars');
                
            else
                cifObj.histVars = {};
                cifObj.histCoeffVars = {};
                histCoeffsVarsTrans = {};
            end
            
            
            
            % Define the functional form of the Conditonal Intensity
            % Function based on how the data was fit.
            cifObj.fitType = fitType;            
            if(isempty(cifObj.histVars))
                if(strcmp(fitType,'poisson'))
                    cifObj.lambdaDelta = simplify(exp(beta*cifObj.varIn)); 
                    cifObj.lambdaDeltaFunction = matlabFunction(cifObj.lambdaDelta,'vars',cifObj.varIn);
                elseif(strcmp(fitType,'binomial'))
                    cifObj.lambdaDelta = simplify(exp(beta*cifObj.varIn)./(1+exp(beta*cifObj.varIn)));
                    cifObj.lambdaDeltaFunction = matlabFunction(cifObj.lambdaDelta,'vars',symvar(cifObj.varIn));
                end
            else
                if(strcmp(fitType,'poisson'))
                    cifObj.lambdaDelta = simplify(exp(beta*cifObj.varIn  + cifObj.histCoeffs*cifObj.histVars)); 
                    cifObj.lambdaDeltaGamma = simplify(exp(beta*cifObj.varIn  + cifObj.histCoeffVars*cifObj.histVars)); 
                    cifObj.lambdaDeltaFunction = matlabFunction(cifObj.lambdaDelta,'vars',[cifObj.varIn; cifObj.histVars]);
                    cifObj.lambdaDeltaGammaFunction = matlabFunction(cifObj.lambdaDeltaGamma,'vars',[cifObj.varIn; cifObj.histVars; histCoeffsVarsTrans]);
                    
                elseif(strcmp(fitType,'binomial'))
                    cifObj.lambdaDelta = simplify(exp(beta*cifObj.varIn  + cifObj.histCoeffs*cifObj.histVars)./(1+exp(beta*cifObj.varIn  + cifObj.histCoeffs*cifObj.histVars)));
                    cifObj.lambdaDeltaGamma = simplify(exp(beta*cifObj.varIn  + cifObj.histCoeffVars*cifObj.histVars)./(1+exp(beta*cifObj.varIn  + cifObj.histCoeffVars*cifObj.histVars)));
                    cifObj.lambdaDeltaFunction = matlabFunction(cifObj.lambdaDelta,'vars',symvar([cifObj.varIn; cifObj.histVars]));
                    cifObj.lambdaDeltaGammaFunction = matlabFunction(cifObj.lambdaDeltaGamma,'vars',symvar([cifObj.varIn; cifObj.histVars; histCoeffsVarsTrans]));
                end
                
                
            end
                
            % Additional Functions needed for decoding
            % The gradient of log(lambda*delta) and the jacobian of 
            % log(lambda*delta)
            cifObj.gradientLambdaDelta = simplify(jacobian(cifObj.lambdaDelta,cifObj.stimVars));
            cifObj.gradientLogLambdaDelta=simplify(jacobian(log(cifObj.lambdaDelta),cifObj.stimVars));
            cifObj.gradientFunction    = matlabFunction(cifObj.gradientLambdaDelta,'vars',[symvar(cifObj.varIn); cifObj.histVars]);
            cifObj.gradientLogFunction = matlabFunction(cifObj.gradientLogLambdaDelta,'vars',[symvar(cifObj.varIn); cifObj.histVars]);
            
            
            cifObj.jacobianLambdaDelta=simplify(jacobian(cifObj.gradientLambdaDelta,cifObj.stimVars));
            cifObj.jacobianFunction = matlabFunction(cifObj.jacobianLambdaDelta,'vars',[symvar(cifObj.varIn); cifObj.histVars]);
            
            cifObj.jacobianLogLambdaDelta=simplify(jacobian(cifObj.gradientLogLambdaDelta,cifObj.stimVars));
            cifObj.jacobianLogFunction = matlabFunction(cifObj.jacobianLogLambdaDelta,'vars',[symvar(cifObj.varIn); cifObj.histVars]);
            
            
            if(and(~isempty(cifObj.histCoeffs),~isempty(cifObj.history)))
                cifObj.LogLambdaDeltaGamma=simplify(log(cifObj.lambdaDeltaGamma));
                cifObj.LogLambdaDeltaGammaFunction = matlabFunction(cifObj.LogLambdaDeltaGamma,'vars',[symvar(cifObj.varIn); cifObj.histVars;histCoeffsVarsTrans]);

                cifObj.gradientLogLambdaDeltaGamma=simplify(jacobian(log(cifObj.lambdaDeltaGamma),cifObj.histCoeffVars));
                cifObj.gradientLambdaDeltaGamma=simplify(jacobian((cifObj.lambdaDeltaGamma),cifObj.histCoeffVars));
                cifObj.gradientLogFunctionGamma = matlabFunction(cifObj.gradientLogLambdaDeltaGamma,'vars',[symvar(cifObj.varIn); cifObj.histVars;histCoeffsVarsTrans]);
                cifObj.gradientFunctionGamma = matlabFunction(cifObj.gradientLambdaDeltaGamma,'vars',[symvar(cifObj.varIn); cifObj.histVars;histCoeffsVarsTrans]);

                cifObj.jacobianLogLambdaDeltaGamma=simplify(jacobian(cifObj.gradientLogLambdaDeltaGamma,cifObj.histCoeffVars));
                cifObj.jacobianLambdaDeltaGamma=simplify(jacobian(cifObj.gradientLambdaDeltaGamma,cifObj.histCoeffVars));
                cifObj.jacobianLogFunctionGamma = matlabFunction(cifObj.jacobianLogLambdaDeltaGamma,'vars',[symvar(cifObj.varIn); cifObj.histVars;histCoeffsVarsTrans]);
                cifObj.jacobianFunctionGamma = matlabFunction(cifObj.jacobianLambdaDeltaGamma,'vars',[symvar(cifObj.varIn); cifObj.histVars;histCoeffsVarsTrans]);

            else
                cifObj.LogLambdaDeltaGamma=[];
                cifObj.LogLambdaDeltaGammaFunction = [];

                cifObj.gradientLogLambdaDeltaGamma=[];
                cifObj.gradientLambdaDeltaGamma=[];
                cifObj.gradientLogFunctionGamma = [];
                cifObj.gradientFunctionGamma = [];

                cifObj.jacobianLogLambdaDeltaGamma=[];
                cifObj.jacobianLambdaDeltaGamma=[];
                cifObj.jacobianLogFunctionGamma = [];
                cifObj.jacobianFunctionGamma = [];

            end
                cifObj.indepVars = symvar(cifObj.lambdaDelta);
            
            
            % Determine the number of variables and make a default string
            % that will be used to evaluate the above functions
            % This is required since functions defined by using the
            % matlabFunction command do not take vector inputs and so each
            % value needs to be passed separatedly. Defining this string
            % now simplifies how we evaluate these functions
            argstr='';
            
            if(length([symvar(cifObj.varIn); cifObj.histVars])==1)
                argstr = 'val';
            else
                for i=1:(length(symvar(cifObj.varIn))+length(cifObj.histVars))
                    if(i==1)
                        argstr = 'val(1)';
                    else 
                        argstr = strcat(argstr,[',val(' num2str(i) ')']);
                    end
                end
                
            end 

            cifObj.argstr = argstr;
            
            
            argstrVarHist='';
            
            if(length([symvar(cifObj.varIn); cifObj.histVars; histCoeffsVarsTrans])==1)
                argstrVarHist = 'val';
            else
                for i=1:(length(symvar(cifObj.varIn))+length(cifObj.histVars)+length(histCoeffsVarsTrans))
                    if(i==1)
                        argstrVarHist = 'val(1)';
                    else 
                        argstrVarHist = strcat(argstrVarHist,[',val(' num2str(i) ')']);
                    end
                end
                
            end 
            
            cifObj.argstrLDGamma = argstrVarHist;
            
            if(~isempty(cifObj.spikeTrain) && ~isempty(cifObj.history))
               cifObj.historyMat = cifObj.history.computeHistory(cifObj.spikeTrain).dataToMatrix;
            else
                cifObj.historyMat = [];
            end
            
        end
        
        function cifObjNew = CIFCopy(cifObj)
%             pause;
%             cifObjNew = CIF(cifObj.b,cifObj.stimVars,cifObj.stimVars,cifObj.fitType);
            %make a new CIF thats super simple
            cifObjNew = CIF([1],['x'],['x'],cifObj.fitType);
            
            %copy parameters from the old cifObj to the new one
            fnames = fields(cifObj);
            for i=1:length(fnames)
                cifObjNew.(fnames{i}) = cifObj.(fnames{i});
            end
            
        end
        
        function setSpikeTrain(cifObj, spikeTrain)
            cifObj.spikeTrain = spikeTrain.nstCopy;
            if(~isempty(cifObj.history))
               cifObj.historyMat = cifObj.history.computeHistory(cifObj.spikeTrain).dataToMatrix;
            else
               cifObj.historyMat = [];
            end
            
            
        end
        function setHistory(cifObj,histObj)
            %Sets the input history object to be the history object that
            %corresponds to this CIF. 
            % histObj: can be of class History or a vector of doubles to be
            %          used in creating a History object
            if(isa(histObj,'History'))
                cifObj.history = History(histObj.windowTimes);
            elseif(isa(histObj,'double'));
                cifObj.history = History(histObj);
            else
                error('History can only be set by passing in a History Object or a vector of windowTimes');
            end
            
            
        end
        
        
        
        function outVal = evalLambdaDelta(cifObj,stimVal,time_index,nst)
            % outVal = evalLambdaDelta(cifObj,stimVal,nst)
            % scalar value of lambda*delta where lambda is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4)
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
           val = [stimVal;histVal];
           evalString = strcat('outVal = cifObj.lambdaDeltaFunction(',cifObj.argstr,');');
           eval(evalString);
        end
        function outVal = evalGradient(cifObj,stimVal,time_index,nst)
            % outVal = evalGradient(cifObj,stimVal,nst)
            % row vector of the gradient of log(lambda*delta) with respect 
            % to the stimulus variables.  
            % The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4)
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
            val = [stimVal;histVal];
            evalString = strcat('outVal = cifObj.gradientFunction(',cifObj.argstr,');');
            eval(evalString);
            
        end
        
         function outVal = evalGradientLog(cifObj,stimVal,time_index,nst)
            % outVal = evalGradient(cifObj,stimVal,nst)
            % row vector of the gradient of log(lambda*delta) with respect 
            % to the stimulus variables.  
            % The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4)
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
            val = [stimVal;histVal];
            evalString = strcat('outVal = cifObj.gradientLogFunction(',cifObj.argstr,');');
            eval(evalString);
            
         end
        
        
        function outVal = evalJacobian(cifObj,stimVal,time_index,nst)
           
            % outVal = evalJacobian(cifObj,stimVal,nst)
            % matrix vector of the jacobian of log(lambda*delta) with 
            % to the stimulus variables. The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            
             if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4)
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            val = [stimVal;histVal];
            evalString = strcat('outVal = cifObj.jacobianFunction(',cifObj.argstr,');');
            eval(evalString);
        end      
        
        function outVal = evalJacobianLog(cifObj,stimVal,time_index,nst)
           
            % outVal = evalJacobian(cifObj,stimVal,nst)
            % matrix vector of the jacobian of log(lambda*delta) with 
            % to the stimulus variables. The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            
             if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4)
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            val = [stimVal;histVal];
            evalString = strcat('outVal = cifObj.jacobianLogFunction(',cifObj.argstr,');');
            eval(evalString);
        end      
        

        
        %%For history parameters
        
        function outVal = evalLDGamma(cifObj,stimVal,time_index,nst,gamma)
            % outVal = evalLambdaDelta(cifObj,stimVal,nst)
            % scalar value of lambda*delta where lambda is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect

            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
           val = [stimVal;histVal;gamma];
           evalString = strcat('outVal = cifObj.lambdaDeltaGammaFunction(',cifObj.argstrLDGamma,');');
           eval(evalString);
        end
         
        function outVal = evalLogLDGamma(cifObj,stimVal,time_index,nst,gamma)
            % outVal = evalLambdaDelta(cifObj,stimVal,nst)
            % scalar value of lambda*delta where lambda is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect

             if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
           val = [stimVal;histVal;gamma];
           evalString = strcat('outVal = cifObj.LogLambdaDeltaGammaFunction(',cifObj.argstrLDGamma,');');
           eval(evalString);
        end
        
        
        function outVal = evalGradientLDGamma(cifObj,stimVal,time_index,nst,gamma)
            % outVal = evalGradient(cifObj,stimVal,nst)
            % row vector of the gradient of log(lambda*delta) with respect 
            % to the stimulus variables.  
            % The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            
            val = [stimVal;histVal;gamma];
            evalString = strcat('outVal = cifObj.gradientFunctionGamma(',cifObj.argstrLDGamma,');');
            eval(evalString);
            
        end
        function outVal = evalGradientLogLDGamma(cifObj,stimVal,time_index,nst,gamma)
            % outVal = evalGradient(cifObj,stimVal,nst)
            % row vector of the gradient of log(lambda*delta) with respect 
            % to the stimulus variables.  
            % The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect
            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end

            
            val = [stimVal;histVal;gamma];
            evalString = strcat('outVal = cifObj.gradientLogFunctionGamma(',cifObj.argstrLDGamma,');');
            eval(evalString);
            
        end
        
        
        
        function outVal = evalJacobianLogLDGamma(cifObj,stimVal,time_index,nst,gamma)
           
            % outVal = evalJacobian(cifObj,stimVal,nst)
            % matrix vector of the jacobian of log(lambda*delta) with 
            % to the stimulus variables. The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect

             if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            val = [stimVal;histVal;gamma];
            evalString = strcat('outVal = cifObj.jacobianLogFunctionGamma(',cifObj.argstrLDGamma,');');
            eval(evalString);
        end
        function outVal = evalJacobianLDGamma(cifObj,stimVal,time_index,nst,gamma)
           
            % outVal = evalJacobian(cifObj,stimVal,nst)
            % matrix vector of the jacobian of log(lambda*delta) with 
            % to the stimulus variables. The gradient is evaluated at the
            % values in stimVal. If there this CIF has history dependence
            % the nspikeTrain nst is used to compute the history effect

            if(nargin<3)
                time_index=[];
                histVal = [];
            end
            
            if(nargin<4 || isempty(nst))
                if(~isempty(time_index) && ~isempty(cifObj.historyMat))
                    histVal=cifObj.historyMat(time_index,:)';
                end
            else
                if(isa(nst,'nspikeTrain'))
                    if(~isempty(cifObj.history))
                        histData=cifObj.history.computeHistory(nst).dataToMatrix;
                        histVal = histData(end,:)';
                    else
                        histVal = [];
                    end
                else
                    error('Second Input must be of class nspikeTrain');
                end
            end
            val = [stimVal;histVal;gamma];
            evalString = strcat('outVal = cifObj.jacobianFunctionGamma(',cifObj.argstrLDGamma,');');
            eval(evalString);
        end
        
        function ans = isSymBeta(cifObj)
          if(isa(cifObj.b,'sym'))
              ans=1;
          else
              ans=0;
          end
        end
       
    end
    
    methods (Static)
        function spikeTrainColl=simulateCIFByThinningFromLambda(lambda,numRealizations,maxTimeRes)%,histCoeffs,histObj)
            % spikeTrainColl=simulateCIFByThinning(lambda,numRealizations,maxTimeRes)
            % Returns a nstColl with numRealization distinct nspikeTrains
            % corresponding to realizations of the point process specified
            % by the conditional intensity function lambda.
            % 
            % lambda: a SignalObj or Covariate that is the CIF time series.
            % numRealizations: number of realizations to return of the
            %                  point process specified by lambda.
            % maxTimeRes: makes sure that only there is only one spike
            %             occurs within the time maxTimeRes. 
            %
            % Note: Currently assumes no history dependence. Needs to be
            %       modified so that a new lambda is determined at each
            %       time step which includes the current spiking activity
            %       in a given realization.
          
%             if(nargin<5)
%                 histObj = [];
%             end
%             
%             if(nargin<4)
%                 histCoeffs = [];
%             end
            
            if(nargin<3)
                maxTimeRes =[];
            end
            Tmax = lambda.maxTime;
            lambdaBound = max(lambda);
            N=ceil(lambdaBound*(1.5*Tmax)); %Expected number of arrivals in interval 1.5*Tmax
            nst=cell(1,numRealizations);
            
            
            for i=1:numRealizations
                u = rand(1,N); %N samples U(0,1)
                w = -log(u)./(lambdaBound); %Exponential rate lambdaBound
                tSpikes = cumsum(w); %time of the spikes
                tSpikes = tSpikes(tSpikes<=Tmax); %keep only in interval of interest

                % Thinning
                
%                 if(and(~isempty(histObj),~isempty(histCoeffs)))
%                     tempnst = nspikeTrain(tSpikes); tempnst.setMinTime(lambda.minTime);
%                     tempnst.setMaxTime(lambda.maxTime);
%                     tempnst.resample(lambda.sampleRate);
%                     histData = histObj.computeHistory(tempnst);
%                     lambdaHist  = SignalObj(lambda.time,exp(histData.dataToMatrix * histCoeffs)); % Assumes poisson lambda 
%                     lambdaProd  = lambda.*lambdaHist;
%                     lambdaBound = max(lambdaProd);
%                     lambdaRatio = lambdaProd.getValueAt(tSpikes)./lambdaBound;
%                 else
                    lambdaRatio = lambda.getValueAt(tSpikes)./lambdaBound;
                                        
%                 end
                u2 = rand(length(lambdaRatio),1);
                %If lambdaRatio is greater than u2 keep spike otherwise throw
                %away
                
                if(~isempty(lambdaRatio))
                    
                    tSpikesThin  = tSpikes(lambdaRatio>=u2);
                else
                    tSpikesThin =[];
                end
                if(isempty(maxTimeRes))
                    nst{i} = nspikeTrain(tSpikesThin);
                    nst{i}.setName(num2str(1));
                else
                    tSpikesThin = unique(ceil(tSpikesThin./maxTimeRes)*maxTimeRes);
                    nst{i} = nspikeTrain(tSpikesThin);
                    nst{i}.setName(num2str(1));
                end
            
            end
            spikeTrainColl=nstColl(nst);
            spikeTrainColl.setMinTime(lambda.minTime);
            spikeTrainColl.setMaxTime(lambda.maxTime);
        end
        
        function [spikeTrainColl, lambda]=simulateCIFByThinning(mu,hist,stim,ens,inputStimSignal,inputEnsSignal,numRealizations,simType)
            % spikeTrainColl=simulateCIF(mu,hist,stim,inputStimSignal,inputEnsSignal,numRealizations)
            % Returns a nstColl with numRealization different nspikeTrain
            % objects. Each nspikeTrain object is one particular
            % realization of the point process defined the input parameters
            % in the following way:
            % lambda*delta = exp(inputTerms)./(1+exp(inputTerms)
            % where inputTerms = (mu + stim*inputStimSignal +
            % hist*spikeTrain + ens*inputEnsSignal)
            %
            % mu: double the indicates the mean rate of the point process
            % hist: a transfer function (tf) object that is convolved with
            %       process spiking activity to determine the history
            %       effect.
            % stim: a transfer function (tf) object that is convolved with
            %       inputStimSignal to determine the stimulus effect
            %
            % ens : a transfer function (tf) object that is convolved with
            %       the inputEnsSignal to determine the ensemble effect
            %
            % inputStimSignal: a SignalObj specifying the stimulation time
            %                  series.
            % inputEnsSignal: a SignalObj specifying the ensemble activity
            %
            % numRealizations: number of nspikeTrains to return. The
            %                  the conditional intensity function will be
            %                  simulated this number of times to generated
            %                  distinct realizations of the point process.
            % <a href="matlab:web('PPSimExample.html', '-helpbrowser')">Example use of simulateCIF</a> 
            %

            if(nargin<8 || isempty(simType))
                simType='binomial';
            end
            if(nargin<7)
                numRealizations =1;
            end
            Ts=hist.Ts;
            if(1/inputStimSignal.sampleRate == hist.Ts && 1/inputStimSignal.sampleRate ==stim.Ts)
                assignin('base','S',stim);
                assignin('base','H',hist);
                assignin('base','E',ens);
                assignin('base','mu',mu);
                assignin('base','Ts',hist.Ts/100);
                assignin('base','TsInt',hist.Ts);
                if(strcmp(simType,'poisson'))
                    simTypeSelect = 1;
                elseif(strcmp(simType,'binomial'))
                    simTypeSelect = 0;
                else
                    error('simType must be either poisson or binomial');
                end
                assignin('base','simTypeSelect',simTypeSelect);
                
                options = simget;
                lambdaData = zeros(length(inputStimSignal.time),numRealizations);
                t=inputStimSignal.time;
                u=[inputStimSignal.data, inputEnsSignal.data];
                assignin('base','t',t);
                assignin('base','u',u);
%                 options.T
                for i=1:numRealizations
                    
                    simOut = sim('PointProcessSimulationThinning','SimulationMode','normal','AbsTol','1e-5',...
                                'SaveState','on','StateSaveName','xout',...
                                'SaveOutput','on','OutputSaveName','yout',...
                                'SaveTime','on','TimeSaveName','tout',...
                                'StopTime', num2str(inputStimSignal.maxTime),...
                                'StartTime', num2str(inputStimSignal.minTime));
                    simOutVars = simOut.who;
                    yout = simOut.get('yout');
                    tout = yout.time;
%                     [tout,~,yout] = sim('PointProcessSimulationThinning',[inputStimSignal.minTime inputStimSignal.maxTime],options,inputStimSignal.dataToStructure, inputEnsSignal.dataToStructure);
                    spikeTimes = tout(yout.signals(1).values>.5);
                    nst{i} = nspikeTrain(spikeTimes);
                    nst{i}.setName(num2str(1));
                    lambdaData(:,i) = interp1(tout, yout.signals(2).values./Ts,inputStimSignal.time);
                end
                spikeTrainColl=nstColl(nst);
                spikeTrainColl.setMinTime(inputStimSignal.minTime);
                spikeTrainColl.setMaxTime(inputStimSignal.maxTime);
                lambda = Covariate(inputStimSignal.time,lambdaData,'\lambda(t|H_t)','time','s','Hz');
            else
                error('History and Stimulus Transfer functions be discrete and have ''Ts'' equal to 1/inputStimSignal.sampleRate');
            end
            
        end
        
        function [spikeTrainColl, lambda]=simulateCIF(mu,hist,stim,ens,inputStimSignal,inputEnsSignal,numRealizations,simType)
            % spikeTrainColl=simulateCIF(mu,hist,stim,inputStimSignal,inputEnsSignal,numRealizations)
            % Returns a nstColl with numRealization different nspikeTrain
            % objects. Each nspikeTrain object is one particular
            % realization of the point process defined the input parameters
            % in the following way:
            % lambda*delta = exp(inputTerms)./(1+exp(inputTerms)
            % where inputTerms = (mu + stim*inputStimSignal +
            % hist*spikeTrain + ens*inputEnsSignal)
            %
            % mu: double the indicates the mean rate of the point process
            % hist: a transfer function (tf) object that is convolved with
            %       process spiking activity to determine the history
            %       effect.
            % stim: a transfer function (tf) object that is convolved with
            %       inputStimSignal to determine the stimulus effect
            %
            % ens : a transfer function (tf) object that is convolved with
            %       the inputEnsSignal to determine the ensemble effect
            %
            % inputStimSignal: a SignalObj specifying the stimulation time
            %                  series.
            % inputEnsSignal: a SignalObj specifying the ensemble activity
            %
            % numRealizations: number of nspikeTrains to return. The
            %                  the conditional intensity function will be
            %                  simulated this number of times to generated
            %                  distinct realizations of the point process.
            % <a href="matlab:web('PPSimExample.html', '-helpbrowser')">Example use of simulateCIF</a> 
            %

            if(nargin<8 || isempty(simType))
                simType='binomial';
            end
            
            if(nargin<7)
                numRealizations =1;
            end
            if(1/inputStimSignal.sampleRate == hist.Ts && 1/inputStimSignal.sampleRate ==stim.Ts)
                assignin('base','S',stim);
                assignin('base','H',hist);
                assignin('base','E',ens);
                assignin('base','mu',mu);
                assignin('base','Ts',stim.Ts);
                if(strcmp(simType,'poisson'))
                    simTypeSelect = 1;
                elseif(strcmp(simType,'binomial'))
                    simTypeSelect = 0;
                else
                    error('simType must be either poisson or binomial');
                end
                assignin('base','simTypeSelect',simTypeSelect);
                
                options = simget;
                lambdaData = zeros(length(inputStimSignal.time),numRealizations);
                for i=1:numRealizations
                    [tout,~,yout] = sim('PointProcessSimulation',[inputStimSignal.minTime inputStimSignal.maxTime],options,inputStimSignal.dataToStructure, inputEnsSignal.dataToStructure);
                    spikeTimes = tout(yout(:,1)>.5);
                    nst{i} = nspikeTrain(spikeTimes);
                    nst{i}.setName(num2str(1));
%                     lambdaData(:,i) = yout(:,2);
                    lambdaData(:,i) = interp1(tout, yout(:,2),inputStimSignal.time);   
                end
                spikeTrainColl=nstColl(nst);
                spikeTrainColl.setMinTime(inputStimSignal.minTime);
                spikeTrainColl.setMaxTime(inputStimSignal.maxTime);
                lambda = Covariate(inputStimSignal.time,lambdaData,'\lambda(t|H_t)','time','s','Hz');
            else
                error('History and Stimulus Transfer functions be discrete and have ''Ts'' equal to 1/inputStimSignal.sampleRate');
            end
            
        end
    
    end
end
