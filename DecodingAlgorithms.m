classdef DecodingAlgorithms
% DECODINGALGORITHMS A class that contains static functions for 
% decoding the hidden states of linear discrete stochastic systems or 
% hybrid linear discrete stochastic systems subject to gaussian noise. 
% The observations can come from either a gaussian observation model 
% or via a point process observation model.
%
% <a href="matlab: methods('DecodingAlgorithms')">methods</a>
% Reference page in Help browser
% <a href="matlab: doc('DecodingAlgorithms')">doc DecodingAlgorithms</a>


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
    end
    
    methods (Static)
          
        %% Point Process Adaptive Filter 
        %PPDecodeFilter takes an object of class CIF describing the
        %conditional intensity function. This routine is more generic since
        %all of the computations for the PPAF are done symbolically based
        %on the CIF object. However, it also means that this version is
        %must slower than the linear version below.
        function [x_p, W_p, x_u, W_u, x_uT,W_uT,x_pT, W_pT, WConvIter] = PPDecodeFilter(A, Q, Px0, dN,lambdaCIFColl,binwidth,x0,Pi0, yT,PiT,estimateTarget,Wconv)  
            % A can be static or can be a different matrix for each time N
            if(nargin<13||isempty(Wconv))
                Wconv =[];
            end
            [C,N]   = size(dN); % N time samples, C cells

            ns=size(A,1); % number of states


            if(nargin<12 || isempty(estimateTarget))
                estimateTarget=0;
            end

            if(nargin<11 || isempty(PiT))
                if(estimateTarget==1)
                    PiT = zeros(size(Q));
                else
                    PiT = 0*diag(ones(ns,1))*1e-6;
                end
            end
            if(nargin<9 || isempty(Pi0))
                Pi0 = zeros(ns,ns);
            end
            if(nargin<10 || isempty(yT))
                yT=[];
                Amat = A;
                Qmat = Q;
                ft   = zeros(size(Amat,2),N);
                PiT = zeros(size(Q));

            else


                PitT= zeros(ns,ns,N);  % Pi(t,T) in Srinivasan et al. 
                QT  = zeros(ns,ns,N);  % The noise covaraince given target observation (Q_t)
                if(estimateTarget==1)
                    PitT(:,:,N)=Q;   % Pi(T,T)=Pi_T + Q_T, setting PiT=0
                else
                    PitT(:,:,N)=PiT+Q;
                end
                PhitT = zeros(ns,ns,N);% phi(t,T) - transition matrix from time T to t
    %             PhiTt = zeros(ns,ns,N);% phi(T,t) - transition matrix from time t to T
                PhitT(:,:,N) = eye(ns,ns); % phi(T,T) = I
                B = zeros(ns,ns,N);    % See Equation 2.21 in Srinivasan et. al

                for n=N:-1:2
                    invA=eye(size(A))/A;
                    % state transition matrix
                    PhitT(:,:,n-1)= invA*PhitT(:,:,n);
    %                 PhiTt(:,:,n)= A^(N-n);

                    % Equation 2.16 in Srinivasan et al. Note there is a typo in the paper. 
                    % This is the correct expression. The term Q_t-1 does not
                    % need to be mulitplied by phi(t-1,t)

                    PitT(:,:,n-1) = invA*PitT(:,:,n)*invA'+Q;



                    if(n<=N)
                        B(:,:,n) = A-(Q/PitT(:,:,n))*A; %Equation 2.21 in Srinivasan et. al
                        QT(:,:,n) = Q-(Q/PitT(:,:,n))*Q';
                    end
                end
    %             PhiTt(:,:,1)= A^(N-1);
                B(:,:,1) = A-(Q/PitT(:,:,1))*A;
                QT(:,:,1) = Q-(Q/PitT(:,:,1))*Q';
                % See Equations 2.23 through 2.26 in Srinivasan et. al
                if(estimateTarget==1)
    %                 beta = [beta ;zeros(ns,C)];
                    for n=1:N
                       psi = B(:,:,n);
                       if(n==N)
                           gammaMat = eye(ns,ns);
                       else
                           gammaMat = (Q/PitT(:,:,n))*PhitT(:,:,n);
                       end
                       Amat(:,:,n) = [psi,gammaMat;
                                      zeros(ns,ns), eye(ns,ns)];
        %                if(n>1)
        %                 tUnc(:,:,n) = tUnc(:,:,n-1)+PhiTt(:,:,n)*Q*PhiTt(:,:,n)';
        %                else
        %                 tUnc(:,:,n) = PhiTt(:,:,n)*Q*PhiTt(:,:,n)';   
        %                end
                       Qmat(:,:,n) = [QT(:,:,n),   zeros(ns,ns);
                                      zeros(ns,ns) zeros(ns,ns)]; 
                    end
                else

                    Amat = B;
                    Qmat = QT;
                    for n=1:N
                        ft(:,n)   = (Q/PitT(:,:,n))*PhitT(:,:,n)*yT;
                    end

                end

            end

            if(nargin<8 || isempty(x0))
                x0=zeros(size(A,2),1);
            end

            if(nargin<7)
                binwidth = .001; % in seconds
            end

            %% 
            % Return values are
            % x_p: state estimates given the past x_k|k-1
            % W_p: error covariance estimates given the past
            % x_u: state updates given the data - x_k|k
            % W_u: error covariance updates given the data

            [C,N]   = size(dN); % N time samples, C cells

              %% Initialize the PPAF
            x_p     = zeros( size(Amat,2), N+1 );
            x_u     = zeros( size(Amat,2), N );
            W_p    = zeros( size(Amat,2),size(Amat,2), N+1 );
            W_u    = zeros( size(Amat,2),size(Amat,2), N );




            if(~isempty(yT))
                if(det(Pi0)==0) % Assume x0 is known exactly

                else %else
                    invPi0 = pinv(Pi0);
                    invPitT= pinv(PitT(:,:,1));
                    Pi0New = pinv(invPi0+invPitT);
                    Pi0New(isnan(Pi0New))=0;
                    x0New  = Pi0New*(invPi0*x0+invPitT*PhitT(:,:,1)*yT);
                    x0=x0New; Pi0 = Pi0New;
                end
            end
            if(~isempty(yT) && estimateTarget==1)
                    x0= [x0;yT]; %simultaneous estimation of target requires state augmentation

            end


            if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
                x_p(:,1)= Amat(:,:,1)*x0;

            else
                invPitT  = pinv(PitT(:,:,1));
    %             invPhitT = pinv(PhitT(:,:,1));
                invA     = pinv(A);
                invPhi0T = pinv(invA*PhitT(:,:,1));
                ut(:,1) = (Q*invPitT)*PhitT(:,:,1)*(yT-invPhi0T*x0);
                [x_p(:,1), W_p(:,:,1)] = DecodingAlgorithms.PPDecode_predict(x0, Pi0, Amat(:,:,min(size(Amat,3),n)), Qmat(:,:,min(size(Qmat,3))));
                x_p(:,1) = x_p(:,1)+ut(:,1);
                W_p(:,:,1) = W_p(:,:,1) + (Q*invPitT)*A*Pi0*A'*(Q*invPitT)';

    %             x_p(:,1)= Amat(:,:,1)*x0 + ft(:,1);


            end
            if(estimateTarget==1 && ~isempty(yT))
               Pi0New = [Pi0, zeros(ns,ns);
                         zeros(ns,ns)  , zeros(ns,ns)];
               W_p(:,:,1) = Amat(:,:,1)*Pi0New*Amat(:,:,1)'+Qmat(:,:,1);      
            elseif(estimateTarget==0 && isempty(yT))

               W_p(:,:,1) = Amat(:,:,1)*Pi0*Amat(:,:,1)'+Qmat(:,:,1);
            end %Otherwise we computed it above.


            for n=1:N
                [x_u(:,n),   W_u(:,:,n)]   = DecodingAlgorithms.PPDecode_update( x_p(:,n), W_p(:,:,n), dN(:,1:n),lambdaCIFColl, binwidth,n);
    %             [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(A,3),n)), Qmat(:,:,min(size(Qmat,3))));

                if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
                    [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(A,3),n)), Qmat(:,:,min(size(Qmat,3))));
                else
                    %ut= Q_{t}\Pi(t,T)^{-1}\phi(t,T)(y_{T}-phi(T,t-1)x_{t-1}
                    if(n<N)
                        ut(:,n+1) = (Q*pinv(PitT(:,:,n+1)))*PhitT(:,:,n+1)*(yT-pinv(PhitT(:,:,n))*x_u(:,n));
        %                 ut(:,n+1) = ut(:,n+1)*delta;
                        [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(A,3),n)), Qmat(:,:,min(size(Qmat,3))));
                        x_p(:,n+1) = x_p(:,n+1)+ut(:,n+1);
                        W_p(:,:,n+1) = W_p(:,:,n+1) + (Q*pinv(PitT(:,:,n+1)))*A*W_u(:,:,n)*A'*(Q*pinv(PitT(:,:,n+1)))';
                    end
                end
                if(n>1 && isempty(Wconv))
                    diffWun = abs(trace(W_u(:,:,n))-W_u(:,:,n-1));
                    mAbsdiffWun = max(max(diffWun));
                    if(mAbsdiffWun<1e-6)
                        Wconv=W_u(:,:,n);
                        WConvIter = n;
                    else
                        WConvIter=[];
                    end
                    
                end
            

            end
            if(~isempty(yT) && estimateTarget==1)
               %decompose the augmented state space into estimates of the state
               %vector and the target position
               x_uT = x_u(ns+1:2*ns,:);
               W_uT = W_u(ns+1:2*ns,ns+1:2*ns,:);
               x_pT = x_p(ns+1:2*ns,:);
               W_pT = W_p(ns+1:2*ns,ns+1:2*ns,:);

               x_u = x_u(1:ns,:);
               W_u = W_u(1:ns,1:ns,:);
               x_p = x_p(1:ns,:);
               W_p = W_p(1:ns,1:ns,:);

            else
               x_uT = [];
               W_uT = [];
               x_pT = [];
               W_pT = [];

            end


        end
        %PPDecodeFilterLinear takes in a linear representation of the
        %conditional intensity terms. These are the terms that are inside
        %the exponential in a Poisson or Binomial description of the CIF.
        %If such a representation is available, use of this routine is
        %recommended because it is much faster.
        function [x_p, W_p, x_u, W_u, x_uT,W_uT,x_pT, W_pT,WConvIter] = PPDecodeFilterLinear(A, Q, dN,mu,beta,fitType,delta,gamma,windowTimes,x0, Pi0, yT,PiT,estimateTarget,Wconv)
        % [x_p, W_p, x_u, W_u] = PPDecodeFilterLinear(CIFType,A, Q, dN,beta,gamma,x0, xT)
        % Point process adaptive filter with the assumption of linear
        % expresion for the conditional intensity functions (see below). If
        % the terms in the conditional intensity function include
        % polynomial powers of a variable for example, these expressions do
        % not hold. Use the PPDecodeFilter instead since it will compute
        % these expressions symbolically. However, because of the matlab
        % symbolic toolbox, it runs much slower than this version.
        %
        % If a final value for xT is given then the approach of Srinivasan
        % et. al (2006) is used for concurrent estimate of the state and
        % the target. This involves state augmentation of the original
        % state space model. If a final value is not specified then the
        % standard Point Process Adaptive Filter of Eden et al (2004) is
        % used instead.
        % 
        % Assumes in both cases that 
        %   x_t = A*x_{t-1} + w_{t}     w_{t} ~ Normal with zero me and
        %                                       covariance Q
        %
        %
        % Paramerters:
        %  
        % A:        The state transition matrix from the x_{t-1} to x_{t}
        %
        % Q:        The covariance of the process noise w_t
        %
        % dN:       A C x N matrix of ones and zeros corresponding to the
        %           observed spike trains. N is the number of time steps in
        %           my code. C is the number of cells
        %
        % mu:       Cx1 vector of baseline firing rates for each cell. In
        %           the CIF expression in 'fitType' description 
        %           mu_c=mu(c);
        %
        % beta:     nsxC matrix of coefficients for the conditional
        %           intensity function. ns is the number of states in x_t 
        %           In the conditional intesity function description below
        %           beta_c = beta(:,c)';
        %
        % fitType: 'poisson' or 'binomial'. Determines how the beta and
        %           gamma coefficients are used to compute the conditional
        %           intensity function.
        %           For the cth cell:
        %           If poisson: lambda*delta = exp(mu_c+beta_c*x + gamma_c*hist_c)
        %           If binomial: logit(lambda*delta) = mu_c+beta_c*x + gamma_c*hist_c
        %
        % delta:    The number of seconds per time step. This is used to compute
        %           th history effect for each spike train using the input
        %           windowTimes and gamma
        %
        % gamma:    length(windowTimes)-1 x C matrix of the history
        %           coefficients for each window in windowTimes. In the 'fitType'
        %           expression above:
        %           gamma_c = gamma(:,c)';
        %           If gamma is a length(windowTimes)-1x1 vector, then the
        %           same history coefficients are used for each cell.
        %
        % windowTimes: Defines the distinct windows of time (in seconds)
        %           that will be computed for each spike train.
        %
        % xT:       Target Position
        %
        % PiT:      Target Uncertainty
        %
        % estimateTarget: By default (==0), it is assumed that that the 
        %                 initial target information is fixed. Set to 1 in order to 
        %                 simultaneously estimate the target location via 
        %                 state augmentation
        %
        %
        %
        % Code for reaching to final target adapted from:
        % L. Srinivasan, U. T. Eden, A. S. Willsky, and E. N. Brown, 
        % "A state-space analysis for reconstruction of goal-directed
        % movements using neural signals.,"
        % Neural computation, vol. 18, no. 10, pp. 2465?2494, Oct. 2006.
        %
        % Point Process Adaptive Filter from 
        % U. T. Eden, L. M. Frank, R. Barbieri, V. Solo, and E. N. Brown, 
        % "Dynamic analysis of neural encoding by point process adaptive
        % filtering.,"
        % Neural computation, vol. 16, no. 5, pp. 971?998, May. 2004.
        
        if(nargin<15||isempty(Wconv))
            Wconv =[];
        end
        [C,N]   = size(dN); % N time samples, C cells
        ns=size(A,1); % number of states
        
        if(nargin<14 || isempty(estimateTarget))
            estimateTarget=0;
        end
        if(nargin<10 || isempty(x0))
           x0=zeros(ns,1);
           
        end
        if(nargin<9 || isempty(windowTimes))
           windowTimes=[]; 
        end
        if(nargin<8 || isempty(gamma))
            gamma=0;
        end
        if(nargin<7 || isempty(delta))
            delta = .001;
        end
        
        if(nargin<13 || isempty(PiT))
            if(estimateTarget==1)
                PiT = zeros(size(Q));
            else
                PiT = 0*diag(ones(ns,1))*1e-6;
            end
        end
        if(nargin<11 || isempty(Pi0))
            Pi0 = zeros(ns,ns);
        end
        if(nargin<12 || isempty(yT))
            yT=[];
            Amat = A;
            Qmat = Q;
            ft   = zeros(size(Amat,2),N);
            PiT = zeros(size(Q));
            
        else
            
            
            PitT= zeros(ns,ns,N);  % Pi(t,T) in Srinivasan et al. 
            QT  = zeros(ns,ns,N);  % The noise covaraince given target observation (Q_t)
            QN =Q(:,:,min(size(Q,3),N));
            if(estimateTarget==1)
                
                PitT(:,:,N)=QN;   % Pi(T,T)=Pi_T + Q_T, setting PiT=0
            else
                PitT(:,:,N)=PiT+QN;
            end
            PhitT = zeros(ns,ns,N);% phi(t,T) - transition matrix from time T to t
%             PhiTt = zeros(ns,ns,N);% phi(T,t) - transition matrix from time t to T
            PhitT(:,:,N) = eye(ns,ns); % phi(T,T) = I
            B = zeros(ns,ns,N);    % See Equation 2.21 in Srinivasan et. al
            
            for n=N:-1:2
                An =A(:,:,min(size(A,3),n));
                Qn =Q(:,:,min(size(Q,3),n));
                
                invA=pinv(An);
                % state transition matrix
                PhitT(:,:,n-1)= invA*PhitT(:,:,n);
%                 PhiTt(:,:,n)= A^(N-n);

                % Equation 2.16 in Srinivasan et al. Note there is a typo in the paper. 
                % This is the correct expression. The term Q_t-1 does not
                % need to be mulitplied by phi(t-1,t)
                
                PitT(:,:,n-1) = invA*PitT(:,:,n)*invA'+Qn;

             

                if(n<=N)
                    
                    B(:,:,n) = An-(Qn*pinv(PitT(:,:,n)))*An; %Equation 2.21 in Srinivasan et. al
                    QT(:,:,n) = Qn-(Qn*pinv(PitT(:,:,n)))*Qn';
                end
            end
            A1=A(:,:,min(size(A,3),1));
            Q1=Q(:,:,min(size(Q,3),1));
            B(:,:,1) = A1-(Q1*pinv(PitT(:,:,1)))*A1;
            QT(:,:,1) = Q1-(Q1*pinv(PitT(:,:,1)))*Q1';

            % See Equations 2.23 through 2.26 in Srinivasan et. al
            if(estimateTarget==1)
                beta = [beta ;zeros(ns,C)];
                for n=1:N
                    An =A(:,:,min(size(A,3),n));
                    Qn =Q(:,:,min(size(Q,3),n));
                    psi = B(:,:,n);
                    if(n==N)
                       gammaMat = eye(ns,ns);
                    else
                       gammaMat = (Qn*pinv(PitT(:,:,n)))*PhitT(:,:,n);
                    end
                    Amat(:,:,n) = [psi,gammaMat;
                                  zeros(ns,ns), eye(ns,ns)];
                    Qmat(:,:,n) = [QT(:,:,n),   zeros(ns,ns);
                                  zeros(ns,ns) zeros(ns,ns)]; 
                end
            else
                
                Amat = B;
                Qmat = QT;
                for n=1:N
                    An =A(:,:,min(size(A,3),n));
                    Qn =Q(:,:,min(size(Q,3),n));
                    ft(:,n)   = (Qn*pinv(PitT(:,:,n)))*PhitT(:,:,n)*yT;
                end

            end

        end
         
        
        minTime=0;
        maxTime=(size(dN,2)-1)*delta;
        
        C=size(dN,1);
        if(~isempty(windowTimes))
            histObj = History(windowTimes,minTime,maxTime);
            HkAll = zeros(size(dN,2),length(windowTimes)-1,C);
            for c=1:C
                nst{c} = nspikeTrain( (find(dN(c,:)==1)-1)*delta);
                nst{c}.setMinTime(minTime);
                nst{c}.setMaxTime(maxTime);
                nst{c}=nst{c}.resample(1/delta);
                HkAll(:,:,c) = histObj.computeHistory(nst{c}).dataToMatrix;
%                 HkAll{c} = histObj.computeHistory(nst{c}).dataToMatrix;
            end
            if(size(gamma,2)==1 && C>1) % if more than 1 cell but only 1 gamma
                gammaNew(:,c) = gamma;
            else
                gammaNew=gamma;
            end
            gamma = gammaNew;
                
        else
            for c=1:C
%                 HkAll{c} = zeros(N,1);
                HkAll(:,:,c) = zeros(N,1);
                gammaNew(c)=0;
            end
            gamma=gammaNew;
            
        end
        if(size(gamma,2)~=C)
            gamma=gamma';
        end
        

        
        %% Initialize the PPAF
        x_p     = zeros( size(Amat,2), N+1 );
        x_u     = zeros( size(Amat,2), N );
        W_p    = zeros( size(Amat,2),size(Amat,2), N+1 );
        W_u    = zeros( size(Amat,2),size(Amat,2), N );
        
        


        if(~isempty(yT))
            if(det(Pi0)==0) % Assume x0 is known exactly
                
            else %else
                invPi0 = pinv(Pi0);
                invPitT= pinv(PitT(:,:,1));
                Pi0New = pinv(invPi0+invPitT);
                Pi0New(isnan(Pi0New))=0;
                x0New  = Pi0New*(invPi0*x0+invPitT*PhitT(:,:,1)*yT);
                x0=x0New; Pi0 = Pi0New;
            end
        end
        if(~isempty(yT) && estimateTarget==1)
                x0= [x0;yT]; %simultaneous estimation of target requires state augmentation
            
        end
        
        
        if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
            x_p(:,1)= Amat(:,:,1)*x0;
           
        else
            invPitT  = pinv(PitT(:,:,1));
%             invPhitT = pinv(PhitT(:,:,1));
            A1 = A(:,:,min(size(A,3),1));
            Q1 = Q(:,:,min(size(Q,3),1));
            invA     = pinv(A1);
            invPhi0T = pinv(invA*PhitT(:,:,1));
            ut(:,1) = (Q1*invPitT)*PhitT(:,:,1)*(yT-invPhi0T*x0);
            [x_p(:,1), W_p(:,:,1)] = DecodingAlgorithms.PPDecode_predict(x0, Pi0, Amat(:,:,min(size(Amat,3),n)), Qmat(:,:,min(size(Qmat,3))));
            x_p(:,1) = x_p(:,1)+ut(:,1);
            W_p(:,:,1) = W_p(:,:,1) + (Q1*invPitT)*A1*Pi0*A1'*(Q1*invPitT)';
                    
%             x_p(:,1)= Amat(:,:,1)*x0 + ft(:,1);

            
        end
        if(estimateTarget==1 && ~isempty(yT))
           Pi0New = [Pi0, zeros(ns,ns);
                     zeros(ns,ns)  , zeros(ns,ns)];
           W_p(:,:,1) = Amat(:,:,1)*Pi0New*Amat(:,:,1)'+Qmat(:,:,1);      
        elseif(estimateTarget==0 && isempty(yT))
            
           W_p(:,:,1) = Amat(:,:,1)*Pi0*Amat(:,:,1)'+Qmat(:,:,1);
        end %Otherwise we computed it above.
        
        HkPerm = permute(HkAll,[2 3 1]);
        clear t;
        for n=1:N


            [x_u(:,n), W_u(:,:,n)] = DecodingAlgorithms.PPDecode_updateLinear(x_p(:,n), W_p(:,:,n), dN,mu,beta,fitType,gamma,HkPerm,n,Wconv);
            % The prediction step is identical to the symbolic implementation since
            % it is independent of the CIF

            if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
                [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(Amat,3),n)), Qmat(:,:,min(size(Qmat,3))),Wconv);
            else
                %ut= Q_{t}\Pi(t,T)^{-1}\phi(t,T)(y_{T}-phi(T,t-1)x_{t-1}
                if(n<N)
                    An = A(:,:,min(size(A,3),n));
                    Qn = Q(:,:,min(size(Q,3),n));
                    invPitT  = pinv(PitT(:,:,n+1));
                    invPhitm1T = pinv(PhitT(:,:,n));
                    ut(:,n+1) = (Qn*invPitT)*PhitT(:,:,n+1)*(yT-invPhitm1T*x_u(:,n));
    %                 ut(:,n+1) = ut(:,n+1)*delta;
                    [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(A,3),n)), Qmat(:,:,min(size(Qmat,3))));
                    x_p(:,n+1) = x_p(:,n+1)+ut(:,n+1);
                    W_p(:,:,n+1) = W_p(:,:,n+1) + (Qn*invPitT)*An*W_u(:,:,n)*An'*(Qn*invPitT)';
                end
            end
%             t(n) = trace(W_u(:,:,n));
%             nSmooth = 100;
%             if(n>nSmooth*3 && isempty(Wconv))
%                 tSig = SignalObj(0:delta:(n-1)*delta, t);
%                 Avals = 1;
%                 Bvals = ones(1,100)./100;
%                 tSig = tSig.filtfilt(Bvals,Avals);
%                 dtSig = tSig.derivative;
%                 diffWun = abs(dtSig.data(n));
%                 if(diffWun<1e-8)
%                     Wconv=W_u(:,:,n);
%                     WConvIter = n;
%                 else
%                     WConvIter=[];
%                 end
%             end
            
%             if(n>1 && isempty(Wconv))
%                 diffWun = abs(trace(W_u(:,:,n))-trace(W_u(:,:,n-1)));
%                 mAbsdiffWun = max(max(diffWun));
% %                 [U,S,V] = svd(diffWun);
% %                 mAbsdiffWun =max(diag(S));
%                 if(mAbsdiffWun<1e-5)
%                     Wconv=W_u(:,:,n);
%                     WConvIter = n;
%                 else
%                     WConvIter=[];
%                 end
%                     
%             end
            WConvIter=[];
        end
%         close all; clear t;
%          figure(10); 
%          for i=1:N
%             t(i)=(trace((W_u(:,:,i)))); 
%          end
%          plot(t); pause
        

        if(~isempty(yT) && estimateTarget==1)
           %decompose the augmented state space into estimates of the state
           %vector and the target position
           x_uT = x_u(ns+1:2*ns,:);
           W_uT = W_u(ns+1:2*ns,ns+1:2*ns,:);
           x_pT = x_p(ns+1:2*ns,:);
           W_pT = W_p(ns+1:2*ns,ns+1:2*ns,:);

           x_u = x_u(1:ns,:);
           W_u = W_u(1:ns,1:ns,:);
           x_p = x_p(1:ns,:);
           W_p = W_p(1:ns,1:ns,:);

        else
           x_uT = [];
           W_uT = [];
           x_pT = [];
           W_pT = [];

        end
        end
      
             %% Point Process Fixed-Interval Smoother
        function  [x_uLag, W_uLag] = PP_fixedIntervalSmoother(A, Q, dN, lags, mu,beta,fitType,delta,gamma,windowTimes,x0, Pi0, yT,PiT,estimateTarget)
        % 
        % Assumes in both cases that 
        %   x_t = A*x_{t-1} + w_{t}     w_{t} ~ Normal with zero me and
        %                                       covariance Q
        %
        %
        % Paramerters:
        %  
        % A:        The state transition matrix from the x_{t-1} to x_{t}
        %
        % Q:        The covariance of the process noise w_t
        %
        % dN:       A C x N matrix of ones and zeros corresponding to the
        %           observed spike trains. N is the number of time steps in
        %           my code. C is the number of cells
        %
        % mu:       Cx1 vector of baseline firing rates for each cell. In
        %           the CIF expression in 'fitType' description 
        %           mu_c=mu(c);
        %
        % beta:     nsxC matrix of coefficients for the conditional
        %           intensity function. ns is the number of states in x_t 
        %           In the conditional intesity function description below
        %           beta_c = beta(:,c)';
        %
        % fitType: 'poisson' or 'binomial'. Determines how the beta and
        %           gamma coefficients are used to compute the conditional
        %           intensity function.
        %           For the cth cell:
        %           If poisson: lambda*delta = exp(mu_c+beta_c*x + gamma_c*hist_c)
        %           If binomial: logit(lambda*delta) = mu_c+beta_c*x + gamma_c*hist_c
        %
        % delta:    The number of seconds per time step. This is used to compute
        %           th history effect for each spike train using the input
        %           windowTimes and gamma
        %
        % gamma:    length(windowTimes)-1 x C matrix of the history
        %           coefficients for each window in windowTimes. In the 'fitType'
        %           expression above:
        %           gamma_c = gamma(:,c)';
        %           If gamma is a length(windowTimes)-1x1 vector, then the
        %           same history coefficients are used for each cell.
        %
        % windowTimes: Defines the distinct windows of time (in seconds)
        %           that will be computed for each spike train.
        %
        % xT:       Target Position
        %
        % PiT:      Target Uncertainty
        %
        % estimateTarget: By default (==0), it is assumed that that the 
        %                 initial target information is fixed. Set to 1 in order to 
        %                 simultaneously estimate the target location via 
        %                 state augmentation
        %
        %
        %
        % Code for reaching to final target adapted from:
        % L. Srinivasan, U. T. Eden, A. S. Willsky, and E. N. Brown, 
        % "A state-space analysis for reconstruction of goal-directed
        % movements using neural signals.,"
        % Neural computation, vol. 18, no. 10, pp. 2465?2494, Oct. 2006.
        %
        % Point Process Adaptive Filter from 
        % U. T. Eden, L. M. Frank, R. Barbieri, V. Solo, and E. N. Brown, 
        % "Dynamic analysis of neural encoding by point process adaptive
        % filtering.,"
        % Neural computation, vol. 16, no. 5, pp. 971?998, May. 2004.
        

        [C,N]   = size(dN); % N time samples, C cells
        ns=size(A,1); % number of states
        
        if(nargin<15 || isempty(estimateTarget))
            estimateTarget=0;
        end
        if(nargin<11 || isempty(x0))
           x0=zeros(ns,1);
           
        end
        if(nargin<10 || isempty(windowTimes))
           windowTimes=[]; 
        end
        if(nargin<9 || isempty(gamma))
            gamma=0;
        end
        if(nargin<8 || isempty(delta))
            delta = .001;
        end
        
        if(nargin<14 || isempty(PiT))
            if(estimateTarget==1)
                PiT = zeros(size(Q));
            else
                PiT = 0*diag(ones(ns,1))*1e-6;
            end
        end
        if(nargin<12 || isempty(Pi0))
            Pi0 = zeros(ns,ns);
        end
        if(nargin<13 || isempty(yT))
            yT=[];
            Amat = A;
            Qmat = Q;
            ft   = zeros(size(Amat,2),N);
            PiT = zeros(size(Q));
            
        else
            
            
            PitT= zeros(ns,ns,N);  % Pi(t,T) in Srinivasan et al. 
            QT  = zeros(ns,ns,N);  % The noise covaraince given target observation (Q_t)
            QN =Q(:,:,min(size(Q,3),N));
            if(estimateTarget==1)
                
                PitT(:,:,N)=QN;   % Pi(T,T)=Pi_T + Q_T, setting PiT=0
            else
                PitT(:,:,N)=PiT+QN;
            end
            PhitT = zeros(ns,ns,N);% phi(t,T) - transition matrix from time T to t
%             PhiTt = zeros(ns,ns,N);% phi(T,t) - transition matrix from time t to T
            PhitT(:,:,N) = eye(ns,ns); % phi(T,T) = I
            B = zeros(ns,ns,N);    % See Equation 2.21 in Srinivasan et. al
            
            for n=N:-1:2
                An =A(:,:,min(size(A,3),n));
                Qn =Q(:,:,min(size(Q,3),n));
                
                invA=pinv(An);
                % state transition matrix
                PhitT(:,:,n-1)= invA*PhitT(:,:,n);
%                 PhiTt(:,:,n)= A^(N-n);

                % Equation 2.16 in Srinivasan et al. Note there is a typo in the paper. 
                % This is the correct expression. The term Q_t-1 does not
                % need to be mulitplied by phi(t-1,t)
                
                PitT(:,:,n-1) = invA*PitT(:,:,n)*invA'+Qn;

             

                if(n<=N)
                    
                    B(:,:,n) = An-(Qn*pinv(PitT(:,:,n)))*An; %Equation 2.21 in Srinivasan et. al
                    QT(:,:,n) = Qn-(Qn*pinv(PitT(:,:,n)))*Qn';
                end
            end
            A1=A(:,:,min(size(A,3),1));
            Q1=Q(:,:,min(size(Q,3),1));
            B(:,:,1) = A1-(Q1*pinv(PitT(:,:,1)))*A1;
            QT(:,:,1) = Q1-(Q1*pinv(PitT(:,:,1)))*Q1';

            % See Equations 2.23 through 2.26 in Srinivasan et. al
            if(estimateTarget==1)
                beta = [beta ;zeros(ns,C)];
                for n=1:N
                    An =A(:,:,min(size(A,3),n));
                    Qn =Q(:,:,min(size(Q,3),n));
                    psi = B(:,:,n);
                    if(n==N)
                       gammaMat = eye(ns,ns);
                    else
                       gammaMat = (Qn*pinv(PitT(:,:,n)))*PhitT(:,:,n);
                    end
                    Amat(:,:,n) = [psi,gammaMat;
                                  zeros(ns,ns), eye(ns,ns)];
                    Qmat(:,:,n) = [QT(:,:,n),   zeros(ns,ns);
                                  zeros(ns,ns) zeros(ns,ns)]; 
                end
            else
                
                Amat = B;
                Qmat = QT;
                for n=1:N
                    An =A(:,:,min(size(A,3),n));
                    Qn =Q(:,:,min(size(Q,3),n));
                    ft(:,n)   = (Qn*pinv(PitT(:,:,n)))*PhitT(:,:,n)*yT;
                end

            end

        end
         
        
        minTime=0;
        maxTime=(size(dN,2)-1)*delta;
        
        C=size(dN,1);
        if(~isempty(windowTimes))
            histObj = History(windowTimes,minTime,maxTime);
            HkAll = zeros(size(dN,2),length(windowTimes)-1,C);
            for c=1:C
                nst{c} = nspikeTrain( (find(dN(c,:)==1)-1)*delta);
                nst{c}.setMinTime(minTime);
                nst{c}.setMaxTime(maxTime);
                nst{c}=nst{c}.resample(1/delta);
                HkAll(:,:,c) = histObj.computeHistory(nst{c}).dataToMatrix;
%                 HkAll{c} = histObj.computeHistory(nst{c}).dataToMatrix;
            end
            if(size(gamma,2)==1 && C>1) % if more than 1 cell but only 1 gamma
                gammaNew(:,c) = gamma;
            else
                gammaNew=gamma;
            end
            gamma = gammaNew;
                
        else
            for c=1:C
%                 HkAll{c} = zeros(N,1);
                HkAll(:,:,c) = zeros(N,1);
                gammaNew(c)=0;
            end
            gamma=gammaNew;
            
        end
        if(size(gamma,2)~=C)
            gamma=gamma';
        end
        

        
        %% Initialize the PPAF
        x_p     = zeros( size(Amat,2), N+1 );
        x_u     = zeros( size(Amat,2), N );
        W_p    = zeros( size(Amat,2),size(Amat,2), N+1 );
        W_u    = zeros( size(Amat,2),size(Amat,2), N );
        
        


        if(~isempty(yT))
            if(det(Pi0)==0) % Assume x0 is known exactly
                
            else %else
                invPi0 = pinv(Pi0);
                invPitT= pinv(PitT(:,:,1));
                Pi0New = pinv(invPi0+invPitT);
                Pi0New(isnan(Pi0New))=0;
                x0New  = Pi0New*(invPi0*x0+invPitT*PhitT(:,:,1)*yT);
                x0=x0New; Pi0 = Pi0New;
            end
        end
        if(~isempty(yT) && estimateTarget==1)
                x0= [x0;yT]; %simultaneous estimation of target requires state augmentation
            
        end
        
        
        if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
            x_p(:,1)= Amat(:,:,1)*x0;
           
        else
            invPitT  = pinv(PitT(:,:,1));
%             invPhitT = pinv(PhitT(:,:,1));
            A1 = A(:,:,min(size(A,3),1));
            Q1 = Q(:,:,min(size(Q,3),1));
            invA     = pinv(A1);
            invPhi0T = pinv(invA*PhitT(:,:,1));
            ut(:,1) = (Q1*invPitT)*PhitT(:,:,1)*(yT-invPhi0T*x0);
            [x_p(:,1), W_p(:,:,1)] = DecodingAlgorithms.PPDecode_predict(x0, Pi0, Amat(:,:,min(size(Amat,3),n)), Qmat(:,:,min(size(Qmat,3))));
            x_p(:,1) = x_p(:,1)+ut(:,1);
            W_p(:,:,1) = W_p(:,:,1) + (Q1*invPitT)*A1*Pi0*A1'*(Q1*invPitT)';
                    
%             x_p(:,1)= Amat(:,:,1)*x0 + ft(:,1);

            
        end
        if(estimateTarget==1 && ~isempty(yT))
           Pi0New = [Pi0, zeros(ns,ns);
                     zeros(ns,ns)  , zeros(ns,ns)];
           W_p(:,:,1) = Amat(:,:,1)*Pi0New*Amat(:,:,1)'+Qmat(:,:,1);      
        elseif(estimateTarget==0 && isempty(yT))
            
           W_p(:,:,1) = Amat(:,:,1)*Pi0*Amat(:,:,1)'+Qmat(:,:,1);
        end %Otherwise we computed it above.
        
        HkPerm = permute(HkAll,[2 3 1]);
        x_uLag = zeros(size(x_u));
        W_uLag = zeros(size(W_u));
        for n=1:N


            [x_u(:,n), W_u(:,:,n)] = DecodingAlgorithms.PPDecode_updateLinear(x_p(:,n), W_p(:,:,n), dN,mu,beta,fitType,gamma,HkPerm,n);
            % The prediction step is identical to the symbolic implementation since
            % it is independent of the CIF
            

            if((estimateTarget==1 && ~isempty(yT)) || isempty(yT))
                [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(Amat,3),n)), Qmat(:,:,min(size(Qmat,3))));
            else
                %ut= Q_{t}\Pi(t,T)^{-1}\phi(t,T)(y_{T}-phi(T,t-1)x_{t-1}
                if(n<N)
                    An = A(:,:,min(size(A,3),n));
                    Qn = Q(:,:,min(size(Q,3),n));
                    invPitT  = pinv(PitT(:,:,n+1));
                    invPhitm1T = pinv(PhitT(:,:,n));
                    ut(:,n+1) = (Qn*invPitT)*PhitT(:,:,n+1)*(yT-invPhitm1T*x_u(:,n));
    %                 ut(:,n+1) = ut(:,n+1)*delta;
                    [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,n), W_u(:,:,n), Amat(:,:,min(size(A,3),n)), Qmat(:,:,min(size(Qmat,3))));
                    x_p(:,n+1) = x_p(:,n+1)+ut(:,n+1);
                    W_p(:,:,n+1) = W_p(:,:,n+1) + (Qn*invPitT)*An*W_u(:,:,n)*An'*(Qn*invPitT)';
                end
            end
            x_K=zeros(ns, lags);
            W_K=zeros(ns,ns,lags);
            LnMinusk=zeros(ns, ns, lags);
            if(n>(lags))
                for k=1:lags
                    
                    LnMinusk(:,:,k)=W_u(:,:,n-k)*A(:,:,min(size(A,3),n-k))'/W_p(:,:,n+1-k);
                    if(k==1)
                        x_K(:,k) = x_u(:,n-k)+LnMinusk(:,:,k)*(x_u(:,n+1-k)-x_p(:,n+1-k));
                        W_K(:,:,k)=W_u(:,:,n-k)+LnMinusk(:,:,k)*(W_u(:,:,n+1-k)-W_p(:,:,n+1-k))*LnMinusk(:,:,k)';
                    else                
                        x_K(:,k) = x_u(:,n-k)+LnMinusk(:,:,k)*(x_K(:,k-1)-x_p(:,n+1-k));
                        W_K(:,:,k)=W_u(:,:,n-k)+LnMinusk(:,:,k)*(W_K(:,:,k-1)-W_p(:,:,n+1-k))*LnMinusk(:,:,k)';
                    end
                    W_K(:,:,k) = 0.5*(W_K(:,:,k)+W_K(:,:,k)');
                end

            end
            x_uLag(:,n)=x_K(:,lags);
            W_uLag(:,:,n)=W_K(:,:,lags);
                

            
        end

        
        
        end
        % PPAF Prediction Step 
        function [x_p, W_p] = PPDecode_predict(x_u, W_u, A, Q,Wconv)
            if((nargin<5)||isempty(Wconv))
                Wconv=[];
            end
                % The PPDecode prediction step 
                    x_p     = A * x_u;

                    if(isempty(Wconv))
                        W_p    = A * W_u * A' + Q;
                    else
                        W_p    = Wconv;
                    end

                    W_p = 0.5*(W_p+W_p');

        end
        % PPAF Update Step 
        %PPDecode_update takes in an object of class CIF
        function [x_u, W_u,lambdaDeltaMat] = PPDecode_update(x_p, W_p, dN,lambdaIn,binwidth,time_index,WuConv)
                % The PPDecode update step that finds the state estimate based on new
                % data

                %Original Code
                if(nargin<7||isempty(WuConv))
                    WuConv=[];
                end
                   clear lambda; 
                if(isa(lambdaIn,'cell'))
                    lambda = lambdaIn;
                elseif(isa(lambdaIn,'CIF'))
                    lambda{1} = lambdaIn;
                else
                    error('Lambda must be a cell of CIFs or a CIF');
                end

                clear gradientMat lambdaDeltaMat;
                sumValVec=zeros(size(W_p,1),1);
                sumValMat=zeros(size(W_p,2),size(W_p,2));
                lambdaDeltaMat = zeros(length(lambda),1);
                for C=1:length(lambda)

                   if(isempty(lambda{C}.historyMat))
                        spikeTimes =(find(dN(C,:)==1)-1)*binwidth;
                        nst = nspikeTrain(spikeTimes);
                        nst.resample(1/binwidth);
                        lambdaDeltaMat(C,1) = lambda{C}.evalLambdaDelta(x_p,time_index,nst);
                        sumValVec = sumValVec+dN(C,end)*lambda{C}.evalGradientLog(x_p,time_index,nst)'-lambda{C}.evalGradient(x_p,time_index,nst)';
                        sumValMat = sumValMat-dN(C,end)*lambda{C}.evalJacobianLog(x_p,time_index,nst)'+lambda{C}.evalJacobian(x_p,time_index,nst)';
                   else % we already have computed the history effect and can just use it - much faster
                        lambdaDeltaMat(C,1) = lambda{C}.evalLambdaDelta(x_p,time_index); 
                        sumValVec = sumValVec+dN(C,end)*lambda{C}.evalGradientLog(x_p,time_index)'-lambda{C}.evalGradient(x_p,time_index)';
                        sumValMat = sumValMat-dN(C,end)*lambda{C}.evalJacobianLog(x_p,time_index)'+lambda{C}.evalJacobian(x_p,time_index)';
                   end


                end


                % Use pinv so that we do a SVD and ignore the zero singular values
                % Sometimes because of the state space model definition and how information
                % is integrated from distinct CIFs the sumValMat is very sparse. This
                % allows us to prevent inverting singular matrices
                if(isempty(WuConv))
%                     invWp = pinv(W_p);
%                     invWu = invWp + sumValMat;
%                     invWu = 0.5*(invWu+invWu');
%                     Wu = pinv(invWu);
                    I=eye(size(W_p));
                    Wu=W_p*(I-(I+sumValMat*W_p)\(sumValMat*W_p));

                    % Make sure that the update covariate is positive definite.
                    W_u=nearestSPD(Wu);
                    W_u=0.5*(W_u+W_u');
                else
                    W_u=0.5*(WuConv+WuConv');
                end
               x_u     = x_p + W_u*(sumValVec);


        end       
        %PPDecode_updateLinear takes in a linear representation of the CIF
        %(much faster)
        function [x_u, W_u,lambdaDeltaMat] = PPDecode_updateLinear(x_p, W_p, dN,mu,beta,fitType,gamma,HkAll,time_index,WuConv)
            C   = size(dN,1); % N time samples, C cells
            if(nargin<10|| isempty(WuConv))
                WuConv=[];
            end
            if(nargin<9 || isempty(time_index))
                time_index=1;
            end
            if(nargin<8 || isempty(HkAll))
                [C,N]=size(dN);
                numWindows = size(HkAll,2);
                
                HkAll = zeros(N,numWindows,C);
%                 HkAll=cell(C,1);
%                 for c=1:C
%                     HkAll{c}=0;
%                 end
            end
            if(nargin<7 || isempty(gamma))
                gamma=zeros(1,C);
            end
            if(nargin<6 || isempty(fitType))
                fitType = 'poisson';
            end

            
            sumValVec=zeros(size(W_p,1),1);
            sumValMat=zeros(size(W_p,2),size(W_p,2));
            lambdaDeltaMat = zeros(C,1);
            if(numel(gamma)==1 && gamma==0)
                gamma = zeros(size(mu))';
            end
            if(strcmp(fitType,'binomial'))
                %  Histtermperm = permute(HkAll,[2 3 1]); need to send it a
                %  permuted version of HkAll
                Histterm = HkAll(:,:,time_index);
                if(size(Histterm,2)~=size(mu,1))
                    Histterm=Histterm';
                end
                linTerm = mu+beta'*x_p + diag(gamma'*Histterm);
                lambdaDeltaMat = exp(linTerm)./(1+exp(linTerm));
                if(any(isnan(lambdaDeltaMat))||any(isinf(lambdaDeltaMat)))
                    indNan = isnan(lambdaDeltaMat);
                    indInf = isinf(lambdaDeltaMat);
                    lambdaDeltaMat(indNan)=1;
                    lambdaDeltaMat(indInf)=1;
                end
                sumValVec=sum(repmat(((dN(:,time_index)-lambdaDeltaMat(:,1)).*(1-lambdaDeltaMat(:,1)))',size(beta,1),1).*beta,2);
                tempVec = ((dN(:,time_index)+(1-2*(lambdaDeltaMat(:,1)))).*(1-(lambdaDeltaMat(:,1))).*(lambdaDeltaMat(:,1)))';
%                 tempVec((tempVec<0))=0;
%                 tempVec((tempVec>1))=1;
                sumValMat = (repmat(tempVec,size(beta,1),1).*beta)*beta';
            elseif(strcmp(fitType,'poisson'))
                Histterm = HkAll(:,:,time_index);
  
                if(~any(gamma~=0))
                    Histterm = Histterm';
                end
                if(size(Histterm,2)~=size(mu,1))
                    Histterm=Histterm';
                end
                linTerm = mu+beta'*x_p + diag(gamma'*Histterm);
                lambdaDeltaMat = exp(linTerm);
                if(any(isnan(lambdaDeltaMat))||any(isinf(lambdaDeltaMat)))
                    indNan = isnan(lambdaDeltaMat);
                    indInf = isinf(lambdaDeltaMat);
                    lambdaDeltaMat(indNan)=1;
                    lambdaDeltaMat(indInf)=1;
                end
                sumValVec=sum(repmat(((dN(:,time_index)-lambdaDeltaMat(:,1)))',size(beta,1),1).*beta,2);
                sumValMat = (repmat(lambdaDeltaMat(:,1)',size(beta,1),1).*beta)*beta';
            end
            if(isempty(WuConv))
                    usePInv=0;
                if(usePInv==1)
                    % Use pinv so that we do a SVD and ignore the zero singular values
                    % Sometimes because of the state space model definition and how information
                    % is integrated from distinct CIFs the sumValMat is very sparse. This
                    % allows us to prevent inverting singular matrices
                    I=eye(size(W_p));
                    Wu=W_p*(I-(I+sumValMat*W_p)\(sumValMat*W_p));
                else
                    I=eye(size(W_p));
                    Wu=W_p*(I-(I+sumValMat*W_p)\(sumValMat*W_p));
                end 
               % Make sure that the update covariance is positive definite.
                W_u=nearestSPD(Wu);
                W_u = Wu;
                W_u=0.5*(W_u+W_u');
            else
                W_u = WuConv;
                W_u=0.5*(W_u+W_u');
            end
%             figure(10); subplot(1,3,2);imagesc(W_u); pause(0.005);
            x_u     = x_p + W_u*(sumValVec);


        end


        %% Point Process Hybrid Filter
        function [S_est, X, W, MU_u, X_s, W_s,pNGivenS]= PPHybridFilterLinear(A, Q, p_ij,Mu0, dN,mu,beta,fitType,binwidth,gamma,windowTimes,x0,Pi0, yT,PiT,estimateTarget,MinClassificationError)

            % General-purpose filter design for neural prosthetic devices.
            % Srinivasan L, Eden UT, Mitter SK, Brown EN.
            % J Neurophysiol. 2007 Oct;98(4):2456-75. Epub 2007 May 23.

            [C,N]   = size(dN); % N time samples, C cells
            nmodels = length(A);
            for s=1:nmodels
                ns(s)=size(A{s},1); % number of states
            end
            nsMax = max(ns);
            if(nargin<17 || isempty(MinClassificationError))
                MinClassificationError=0; %0: chooses the most probable discrete state estimate and take the
                                          %   probability weighted average
                                          %   of the continous states. This
                                          %   is the approximate MMSE
                                          %   filter.
                                          %1: takes the most likely discrete state estimate and also the 
                                          %   continuous states
                                          %   corresponding to the most
                                          %   likely discrete state model.
                                          %   This is approximately the
                                          %   Maximum Likelihood Filter
            end
            if(nargin<16 || isempty(estimateTarget))
                estimateTarget=0;
            end

            if(nargin<15 || isempty(PiT))
                for s=1:nmodels
                    if(estimateTarget==1)
                        PiT{s} = zeros(size(Q{s}));
                    else
                        PiT{s} = 0*diag(ones(ns(s),1))*1e-6;
                    end
                end
            end
            if(nargin<13 || isempty(Pi0))
                for s=1:nmodels
                    Pi0{s}(:,:) = zeros(ns(s),ns(s));
                end
            end
            if(nargin<14 || isempty(yT))
                for s=1:nmodels
                    yT{s}=[];
                    Amat{s} = A{s};
                    Qmat{s} = Q{s};
                    ft{s}   = zeros(size(Amat{s},2),N);
                    PiT{s} = zeros(size(Q{s}));
                    betaNew{s} = beta;
                end
                beta = betaNew;
            else
                Pi0new=cell(1,nmodels);
                for s=1:nmodels

                    PitT{s}= zeros(ns(s),ns(s),N);  % Pi(t,T) in Srinivasan et al. 
                    QT{s}  = zeros(ns(s),ns(s),N);  % The noise covaraince given target observation (Q_t)
                    if(estimateTarget==1)
                        PitT{s}(:,:,N)=Q{s};   % Pi(T,T)=Pi_T + Q_T, setting PiT=0
                    else
                        PitT{s}(:,:,N)=PiT{s}+Q{s};
                    end
                    PhitT{s} = zeros(ns(s),ns(s),N);% phi(t,T) - transition matrix from time T to t
        %             PhiTt = zeros(ns,ns,N);% phi(T,t) - transition matrix from time t to T
                    PhitT{s}(:,:,N) = eye(ns(s),ns(s)); % phi(T,T) = I
                    B{s} = zeros(ns(s),ns(s),N);    % See Equation 2.21 in Srinivasan et. al

                    for n=N:-1:2
                        if(rcond(A{s})<1000*eps)
                            invA=pinv(A{s});
                        else
                            invA=eye(size(A{s}))/A{s};
                        end
                        % state transition matrix
                        PhitT{s}(:,:,n-1)= invA*PhitT{s}(:,:,n);
        %                 PhiTt(:,:,n)= A^(N-n);

                        % Equation 2.16 in Srinivasan et al. Note there is a typo in the paper. 
                        % This is the correct expression. The term Q_t-1 does not
                        % need to be mulitplied by phi(t-1,t)

                        PitT{s}(:,:,n-1) = invA*PitT{s}(:,:,n)*invA'+Q{s};


                        if(n<=N)
                            B{s}(:,:,n) = A{s}-(Q{s}*pinv(PitT{s}(:,:,n)))*A{s}; %Equation 2.21 in Srinivasan et. al
                            QT{s}(:,:,n) = Q{s}-(Q{s}*pinv(PitT{s}(:,:,n)))*Q{s}';
                        end
                    end
        %             PhiTt(:,:,1)= A^(N-1);
                    B{s}(:,:,1) = A{s}-(Q{s}*pinv(PitT{s}(:,:,1)))*A{s};
                    QT{s}(:,:,1) = Q{s}-(Q{s}*pinv(PitT{s}(:,:,1)))*Q{s}';
                    betaNew{s} = beta;
                    % See Equations 2.23 through 2.26 in Srinivasan et. al
                    if(estimateTarget==1)

                        betaNew{s} = [beta ;zeros(ns(s),C)];
                        for n=1:N
                           psi = B{s}(:,:,n);
                           if(n==N)
                               gammaMat = eye(ns(s),ns(s));
                           else
                               gammaMat = (Q{s}*pinv(PitT{s}(:,:,n)))*PhitT{s}(:,:,n);
                           end
                           Amat{s}(:,:,n) = [psi,gammaMat;
                                          zeros(ns(s),ns(s)), eye(ns(s),ns(s))];
            %                if(n>1)
            %                 tUnc(:,:,n) = tUnc(:,:,n-1)+PhiTt(:,:,n)*Q*PhiTt(:,:,n)';
            %                else
            %                 tUnc(:,:,n) = PhiTt(:,:,n)*Q*PhiTt(:,:,n)';   
            %                end
                           Qmat{s}(:,:,n) = [QT{s}(:,:,n),   zeros(ns(s),ns(s));
                                          zeros(ns(s),ns(s)) zeros(ns(s),ns(s))]; 



                        end

                       Pi0new{s} = [Pi0{s},  zeros(ns(s),ns(s));
                                    zeros(ns(s),ns(s)) zeros(ns(s),ns(s))]; 

                    else

                        Amat{s} = B{s};
                        Qmat{s} = QT{s};
                        for n=1:N
                            ft{s}(:,n)   = (Q{s}*pinv(PitT{s}(:,:,n)))*PhitT{s}(:,:,n)*yT{s};
                        end

                    end

                end
                if(estimateTarget==1)
                    Pi0 = Pi0new;

                end
                beta = betaNew;
            end

            if(nargin<12 || isempty(x0))
                for s=1:nmodels
                    x0{s}=zeros(size(Amat{s},2),1);
                end
            end

            if(nargin<9)
                binwidth=0.001; %1 msec
            end

            if(isa(A,'cell'))
                dimMat = zeros(1,length(Amat));
                X_u = cell(1,length(Amat));
                W_u = cell(1,length(Amat));
                X_p = cell(1,length(Amat));
                W_p = cell(1,length(Amat));
                ind = cell(1,length(Amat));
                ut = cell(1,length(Amat));

                for i=1:length(Amat)
                    lambdaDeltaMat{i} = zeros(size(dN));
                    X_u{i} = zeros(size(Amat{i},1), size(dN,2));
                    X_p{i} = zeros(size(Amat{i},1), size(dN,2)+1);
                    W_u{i} = zeros(size(Amat{i},1), size(Amat{i},1), size(dN,2));
                    W_p{i} = zeros(size(Amat{i},1), size(Amat{i},1), size(dN,2)+1);
                    dimMat(i) = size(Amat{i},2);
                    W_u{i}(:,:,1) =Pi0{i};
                    ind{i} = 1:dimMat(i);
                    ut{i} = zeros(size(Amat{i},1), size(dN,2));
                end
            end

            maxDim = max(dimMat);
    %         nmodels = length(Amat);
    %         lambdaCIFColl = CIFColl(lambda);

            minTime=0;
            maxTime=(size(dN,2)-1)*binwidth;

            C=size(dN,1);

            if(nargin<11 || isempty(windowTimes))
                 for c=1:C
                    HkAll(:,:,c) = zeros(N,1);
                    gammaNew(c)=0;
                end
                gamma=gammaNew;
            else
                histObj = History(windowTimes,minTime,maxTime);
                for c=1:C
                    nst{c} = nspikeTrain( (find(dN(c,:)==1)-1)*binwidth);
                    nst{c}.setMinTime(minTime);
                    nst{c}.setMaxTime(maxTime);
                    nst{c}=nst{c}.resample(1/delta);
                    HkAll(:,:,c) = histObj.computeHistory(nst{c}).dataToMatrix;
                end
                if(size(gamma,2)==1 && C>1) % if more than 1 cell but only 1 gamma
                    gammaNew(:,c) = gamma;
                end
                gamma = gammaNew;

            end
            % Overall estimates of Hybrid filter
            X = zeros(maxDim, size(dN,2));         % Estimated Trajectories
            W = zeros(maxDim, maxDim, size(dN,2)); % Covariance of estimate
            % Individual Model Estimates
            for i=1:nmodels
                X_s{i} = X;    % Individual Model Estimates
                W_s{i} = W;    % Individual Model Covariances
            end
            % Model probabilities 
            MU_u = zeros(nmodels,size(dN,2));   % P(s_k | H_k+1) % updated state probabilities
            MU_p  = zeros(nmodels,size(dN,2));   % P(s_k | H_k)  % prediction state probabilities
            pNGivenS = zeros(nmodels,size(dN,2));


            %mu_0|1 = mu_0|0;
            if(isempty(Mu0))
                Mu0 = ones(nmodels,1)*1/nmodels;
            elseif(size(Mu0,1)==nmodels && size(Mu0,2)==1)
                Mu0 = Mu0;
            elseif(size(Mu0,1)==1 && size(Mu0,2)==nmodels)
                Mu0 = Mu0';
            else
                error('Mu0 must be a column or row vector with the same number of dimensions as the number of states');
            end
            for s=1:nmodels
                [X_p{s}(ind{s},1),W_p{s}(ind{s},ind{s},1)] = DecodingAlgorithms.PPDecode_predict(x0{s}(ind{s}), Pi0{s}(ind{s},ind{s}), Amat{s}(ind{s},ind{s},min(size(Amat{s},3),1)), Qmat{s}(:,:,min(size(Qmat{s},3))));
              
                if((estimateTarget==0 && ~isempty(yT{s})))               
                    invA= pinv(Amat{s}(:,:,min(size(Amat,3),1)));
                    ut{s}(:,1) = (Q{s}*pinv(PitT{s}(:,:,1)))*PhitT{s}(:,:,1)*(yT{s}-pinv(invA*PhitT{s}(:,:,1))*x0{s});
                    X_p{s}(ind{s},1) = X_p{s}(ind{s},1)+ut{s}(:,1);
                    W_p{s}(ind{s},ind{s},1) =W_p{s}(ind{s},ind{s},1) + (Q{s}*pinv(PitT{s}(:,:,1)))*A{s}*Pi0{s}*A{s}'*(Q{s}*pinv(PitT{s}(:,:,1)))';
                end
            end

    %            [~, S_est(1)] = max(MU_p(:,1)); %Most likely current state

                %State transition Probabilities must integrate to 1
                sumVal = sum(p_ij,2);
                if(any(sumVal~=1))
                    error('State Transition probability matrix must sum to 1 along each row');
                end
            %% 9 Steps
            % Filtering steps.
            HkPerm=permute(HkAll,[2 3 1]);
            for k = 1:(size(dN,2))

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 1 - p(s_k | H_k) = Sum(p(s_k|s_k-1)*p(s_k-1|H_k))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %MU(:,k)_p= [ p(s_k=1 | H_k)] is a vector of the probabilities
                %           [ p(s_k=2 | H_k)]
                %               .
                %               .
                %               .
                %           [ p(s_k=N | H_k)]
                % thus it is an prediction of the discrete state at time k given all
                % of the neural firing up to k-1 as summarized in H_k
                %
                % Whereas 
                % MU_u(:,k)=[ p(s_k=1 | H_k+1)] is a vector of the probabilities
                %           [ p(s_k=2 | H_k+1)]
                %               .
                %               .
                %               .
                %           [ p(s_k=N | H_k+1)]
                % The s suffix indicates that this is a "smoothed" estimate of
                % the state given the firing up to time k summarized in H_k+1
                if k==1
                    MU_p(:,k) = p_ij'*Mu0;          %state probability prediction equation
                else
                    MU_p(:,k) = p_ij'*MU_u(:,k-1);  
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 2 - p(s_k-1 | s_k, H_k)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % This is a matrix with i,j entry indicating the probability that 
                % s_k-1 = j given than s_k = i
                %
                % MU_p is the normalization factor. The first column of the
                % matrix of probabilities is:
                %
                % P(s_k-1=1 | s_k=1,H_k) ~ P(s_k=1|s_k-1=1,H_k)*P(s_k-1=1|H_k)
                % P(s_k-1=1 | s_k=2,H_k) ~ P(s_k=2|s_k-1=1,H_k)*P(s_k-1=1|H_k)
                % 
                % And the second columns ... etc
                %
                % 
                % P(s_k-1=2 | s_k=1,H_k) ~ P(s_k=1|s_k-1=2,H_k)*P(s_k-1=1|H_k)
                % P(s_k-1=2 | s_k=2,H_k) ~ P(s_k=2|s_k-1=2,H_k)*P(s_k-1=1|H_k)

                if(k==1)
                    p_ij_s= p_ij.*(Mu0*ones(1,nmodels));%.*(ones(nmodels,1)*(1./MU_p(:,k))');
                else
                    p_ij_s= p_ij.*(MU_u(:,k-1)*ones(1,nmodels));%.*(ones(nmodels,1)*(1./MU_p(:,k))');
                end
    %          
                 % To avoid any numerical issues with roundoff, we normalize to
                 % 1 again
                 normFact = repmat(sum(p_ij_s,1),[nmodels 1]); %Every column must sum to 1

                 p_ij_s = p_ij_s./normFact;
    %              for i=1:length(normFact)
    %                  if(normFact(i)~=0)
    %                     p_ij_s(:,i) = p_ij_s(:,i)./normFact(i);
    %                  else %reset all the states to be equally likely (each row must sum to 1)
    %                     p_ij_s(:,i) = 1/nmodels*ones(nmodels,1);
    %                  end
    %              end         
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 3 - Approximate p(x_k-1 | s_k, H_k) with Gaussian 
                % approximation to Mixtures of Gaussians
                % Calculate the mixed state mean for each filter
                % This will be the initial states for the update step of the
                % Point Process Adaptive Filter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    for j = 1:nmodels
                        for i = 1:nmodels
                            if(k>1)
                                X_s{j}(ind{i},k) = X_s{j}(ind{i},k) + X_u{i}(ind{i},k-1)*p_ij_s(i,j);
                            else
                                X_s{j}(ind{i},k) = X_s{j}(ind{i},k) + x0{i}(ind{i})*p_ij_s(i,j); 
                            end
                        end
                    end

                        % Calculate the mixed state covariance for each filter

                    for j = 1:nmodels
                        for i = 1:nmodels
                            if(k>1)
                                W_s{j}(ind{i},ind{i},k) = W_s{j}(ind{i},ind{i},k) + (W_u{i}(ind{i},ind{i},k-1) + (X_u{i}(ind{i},k-1)-X_s{j}(ind{i},k))*(X_u{i}(ind{i},k-1)-X_s{j}(ind{i},k))')*p_ij_s(i,j);
                            else
                                W_s{j}(ind{i},ind{i},k) = W_s{j}(ind{i},ind{i},k) + (Pi0{i}(ind{i},ind{i})+ (x0{i}(ind{i})-X_s{j}(ind{i},k))*(x0{i}(ind{i})-X_s{j}(ind{i},k))')*p_ij_s(i,j);
                            end
                        end
                    end

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 4 - Approximate p(x_k+1 |s_k+1,n_k+1,H_k+1)
               % Uses a bank of nmodel point process filters
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            k
    
               for s=1:nmodels

                   % Prediction Step
                   [X_p{s}(ind{s},k),W_p{s}(ind{s},ind{s},k)] = DecodingAlgorithms.PPDecode_predict(X_s{s}(ind{s},k), W_s{s}(ind{s},ind{s},k), Amat{s}(:,:,min(size(Amat,3),k)), Qmat{s}(:,:,min(size(Qmat{s},3))));

                   if(estimateTarget==0 && ~isempty(yT{s}))
                       if(k>1)
                            ut{s}(:,k) = (Q{s}*pinv(PitT{s}(:,:,k)))*PhitT{s}(:,:,k)*(yT{s}-pinv(PhitT{s}(:,:,k-1))*X_s{s}(ind{s},k));
                       else
                           invA = pinv(A{s}(:,:,min(size(A{s},3),1)));
                           ut{s}(:,k) = (Q{s}*pinv(PitT{s}(:,:,1)))*PhitT{s}(:,:,1)*(yT{s}-pinv(invA*PhitT{s}(:,:,1))*X_s{s}(ind{s},k));
                       end
                       X_p{s}(ind{s},k) = X_p{s}(ind{s},k)+ut{s}(:,k);
                       W_p{s}(ind{s},ind{s},k) =W_p{s}(ind{s},ind{s},k) + (Q{s}*pinv(PitT{s}(:,:,k)))*A{s}*W_s{s}(ind{s},ind{s},k)*A{s}'*(Q{s}*pinv(PitT{s}(:,:,k)))';

                    end

                   % Update Step
                   % Fold in the neural firing in the current time step
                   [X_u{s}(ind{s},k),W_u{s}(ind{s},ind{s},k),lambdaDeltaMat{s}(:,k)] = DecodingAlgorithms.PPDecode_updateLinear(X_p{s}(ind{s},k),squeeze(W_p{s}(ind{s},ind{s},k)),dN,mu,beta{s}(ind{s},:),fitType,gamma,HkPerm,k);


               end
    %            pause;
    %            close all; plot(lambdaDeltaMat{1}(:,k),'k.'); hold on; plot(lambdaDeltaMat{2}(:,k),'b*');
    %            
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 5 - p(n_k | s_k, H_k) using Laplace approximation
               % See General-purpose filter design for neural prosthetic devices.
               % Srinivasan L, Eden UT, Mitter SK, Brown EN.
               % J Neurophysiol. 2007 Oct;98(4):2456-75. Epub 2007 May 23.
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


               for s=1:nmodels

                 tempPdf = sqrt(det(W_u{s}(:,:,k)))./sqrt(det(W_p{s}(:,:,k)))*prod(exp(dN(:,k).*log(lambdaDeltaMat{s}(:,k))-lambdaDeltaMat{s}(:,k)));
                 pNGivenS(s,k) = tempPdf;
               end
               tempData = pNGivenS(:,k);
               tempData(isinf(tempData))=0;
               pNGivenS(:,k) = tempData;

               normFact = sum(pNGivenS(:,k));
               if(normFact~=0 && ~isnan(normFact))
                pNGivenS(:,k)=pNGivenS(:,k)./sum(pNGivenS(:,k));
               else

                   if(k>1)
                       pNGivenS(:,k) = pNGivenS(:,k-1);
                   else
                       pNGivenS(:,k) = .5*ones(nmodels,1);
                   end
               end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 6 - Calculate p(s_k | n_k, H_k) = p(s_k | H_k+1)
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               pSGivenN(:,k) = MU_p(:,k).*pNGivenS(:,k);



               %Normalization Factor
               normFact = sum(pSGivenN(:,k));
               if(normFact~=0 && ~isnan(normFact))
                pSGivenN(:,k) = pSGivenN(:,k)./sum(pSGivenN(:,k));
               else
                   if(k>1)
                        pSGivenN(:,k) = pSGivenN(:,k-1);
                   else
                        pSGivenN(:,k) = Mu0;
                   end

               end


               MU_u(:,k) = pSGivenN(:,k); %estimate of s_k given data up to k


               [~, S_est(k)] = max(MU_u(:,k)); %Most likely current state

               if(MinClassificationError==1)

                   s= S_est(k);
                   X(ind{s},k) = X_u{s}(ind{s},k);
                   W(ind{s},ind{s},k) = W_u{s}(ind{s},ind{s},k);

               else %Minimize the mean squared error

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   % Step 7 - Calculate p(x_k | n_k, H_k) - using gaussian
                   % approximation to mixture of gaussians 
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   for s=1:nmodels
                       X(ind{s},k) = X(ind{s},k)+MU_u(s,k)*X_u{s}(ind{s},k); 
                   end
                   for s=1:nmodels
                       W(ind{s},ind{s},k) =  W(ind{s},ind{s},k) +MU_u(s,k)*(W_u{s}(ind{s},ind{s},k) + (X_u{s}(ind{s},k)-X(ind{s},k))*(X_u{s}(ind{s},k)-X(ind{s},k))');
                   end
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               end
            end


            %Solution to the known target problem. Run the Hybrid Filter
            %Forward to determine the most likely states.... this is the
            %sequence of the A dynamics matrices to use in the computation of
            %the Target Reach Model of Srinivasan et al. Compute PiT

         end   
        function [S_est, X, W, MU_s, X_s, W_s,pNGivenS]  = PPHybridFilter(A, Q, p_ij,Mu0,dN,lambdaCIFColl,binwidth,x0,Pi0, yT,PiT,estimateTarget,MinClassificationError)


            % General-purpose filter design for neural prosthetic devices.
            % Srinivasan L, Eden UT, Mitter SK, Brown EN.
            % J Neurophysiol. 2007 Oct;98(4):2456-75. Epub 2007 May 23.

            [C,N]   = size(dN); % N time samples, C cells
            nmodels = length(A);
            for s=1:nmodels
                ns(s)=size(A{s},1); % number of states
            end
            nsMax = max(ns);
            if(nargin<13 || isempty(MinClassificationError))
                MinClassificationError=0; %Minimum mean square error state estimate. By default do maximum likelihood
            end
            if(nargin<12 || isempty(estimateTarget))
                estimateTarget=0;
            end

            if(nargin<11 || isempty(PiT))
                for s=1:nmodels
                    if(estimateTarget==1)
                        PiT{s} = zeros(size(Q{s}));
                    else
                        PiT{s} = 0*diag(ones(ns(s),1))*1e-6;
                    end
                end
            end
            if(nargin<9 || isempty(Pi0))
                for s=1:nmodels
                    Pi0{s}(:,:) = zeros(ns(s),ns(s));
                end
            end
            if(nargin<10 || isempty(yT))
                for s=1:nmodels
                    yT{s}=[];
                    Amat{s} = A{s};
                    Qmat{s} = Q{s};
                    ft{s}   = zeros(size(Amat{s},2),N);
                    PiT{s} = zeros(size(Q{s}));
                end
         
            else
                Pi0new=cell(1,nmodels);
                for s=1:nmodels

                    PitT{s}= zeros(ns(s),ns(s),N);  % Pi(t,T) in Srinivasan et al. 
                    QT{s}  = zeros(ns(s),ns(s),N);  % The noise covaraince given target observation (Q_t)
                    if(estimateTarget==1)
                        PitT{s}(:,:,N)=Q{s};   % Pi(T,T)=Pi_T + Q_T, setting PiT=0
                    else
                        PitT{s}(:,:,N)=PiT{s}+Q{s};
                    end
                    PhitT{s} = zeros(ns(s),ns(s),N);% phi(t,T) - transition matrix from time T to t
        %             PhiTt = zeros(ns,ns,N);% phi(T,t) - transition matrix from time t to T
                    PhitT{s}(:,:,N) = eye(ns(s),ns(s)); % phi(T,T) = I
                    B{s} = zeros(ns(s),ns(s),N);    % See Equation 2.21 in Srinivasan et. al

                    for n=N:-1:2
                        if(rcond(A{s})<1000*eps)
                            invA=pinv(A{s});
                        else
                            invA=eye(size(A{s}))/A{s};
                        end
                        % state transition matrix
                        PhitT{s}(:,:,n-1)= invA*PhitT{s}(:,:,n);
        %                 PhiTt(:,:,n)= A^(N-n);

                        % Equation 2.16 in Srinivasan et al. Note there is a typo in the paper. 
                        % This is the correct expression. The term Q_t-1 does not
                        % need to be mulitplied by phi(t-1,t)

                        PitT{s}(:,:,n-1) = invA*PitT{s}(:,:,n)*invA'+Q{s};


                        if(n<=N)
                            B{s}(:,:,n) = A{s}-(Q{s}*pinv(PitT{s}(:,:,n)))*A{s}; %Equation 2.21 in Srinivasan et. al
                            QT{s}(:,:,n) = Q{s}-(Q{s}*pinv(PitT{s}(:,:,n)))*Q{s}';
                        end
                    end
        %             PhiTt(:,:,1)= A^(N-1);
                    B{s}(:,:,1) = A{s}-(Q{s}*pinv(PitT{s}(:,:,1)))*A{s};
                    QT{s}(:,:,1) = Q{s}-(Q{s}*pinv(PitT{s}(:,:,1)))*Q{s}';
                    % See Equations 2.23 through 2.26 in Srinivasan et. al
                    if(estimateTarget==1)

                   
                        for n=1:N
                           psi = B{s}(:,:,n);
                           if(n==N)
                               gammaMat = eye(ns(s),ns(s));
                           else
                               gammaMat = (Q{s}*pinv(PitT{s}(:,:,n)))*PhitT{s}(:,:,n);
                           end
                           Amat{s}(:,:,n) = [psi,gammaMat;
                                          zeros(ns(s),ns(s)), eye(ns(s),ns(s))];
          
                           Qmat{s}(:,:,n) = [QT{s}(:,:,n),   zeros(ns(s),ns(s));
                                          zeros(ns(s),ns(s)) zeros(ns(s),ns(s))]; 


                        end

                       Pi0new{s} = [Pi0{s},  zeros(ns(s),ns(s));
                                    zeros(ns(s),ns(s)) zeros(ns(s),ns(s))]; 

                    else

                        Amat{s} = B{s};
                        Qmat{s} = QT{s};
                        for n=1:N
                            ft{s}(:,n)   = (Q{s}*pinv(PitT{s}(:,:,n)))*PhitT{s}(:,:,n)*yT{s};
                        end

                    end

                end
                if(estimateTarget==1)
                    Pi0 = Pi0new;

                end
            end

            if(nargin<8 || isempty(x0))
                for s=1:nmodels
                    x0{s}=zeros(size(Amat{s},2),1);
                end
            end

            if(nargin<7)
                binwidth=0.001; %1 msec
            end

            if(isa(A,'cell'))
                dimMat = zeros(1,length(Amat));
                X_u = cell(1,length(Amat));
                W_u = cell(1,length(Amat));
                X_p = cell(1,length(Amat));
                W_p = cell(1,length(Amat));
                ind = cell(1,length(Amat));
                ut = cell(1,length(Amat));

                for i=1:length(Amat)
                    lambdaDeltaMat{i} = zeros(size(dN));
                    X_u{i} = zeros(size(Amat{i},1), size(dN,2));
                    X_p{i} = zeros(size(Amat{i},1), size(dN,2)+1);
                    W_u{i} = zeros(size(Amat{i},1), size(Amat{i},1), size(dN,2));
                    W_p{i} = zeros(size(Amat{i},1), size(Amat{i},1), size(dN,2)+1);
                    dimMat(i) = size(Amat{i},2);
                    W_u{i}(:,:,1) =Pi0{i};
                    ind{i} = 1:dimMat(i);
                    ut{i} = zeros(size(Amat{i},1), size(dN,2));
                end
            end

            maxDim = max(dimMat);
    %         nmodels = length(Amat);
    %         lambdaCIFColl = CIFColl(lambda);

            minTime=0;
            maxTime=(size(dN,2)-1)*binwidth;

            C=size(dN,1);

            % Overall estimates of Hybrid filter
            X = zeros(maxDim, size(dN,2));         % Estimated Trajectories
            W = zeros(maxDim, maxDim, size(dN,2)); % Covariance of estimate
            % Individual Model Estimates
            for i=1:nmodels
                X_s{i} = X;    % Individual Model Estimates
                W_s{i} = W;    % Individual Model Covariances
            end
            % Model probabilities 
            MU_u = zeros(nmodels,size(dN,2));   % P(s_k | H_k+1) % updated state probabilities
            MU_p  = zeros(nmodels,size(dN,2));   % P(s_k | H_k)  % prediction state probabilities
            pNGivenS = zeros(nmodels,size(dN,2));


            %mu_0|1 = mu_0|0;
            if(isempty(Mu0))
                Mu0 = ones(nmodels,1)*1/nmodels;
            elseif(size(Mu0,1)==nmodels && size(Mu0,2)==1)
                Mu0 = Mu0;
            elseif(size(Mu0,1)==1 && size(Mu0,2)==nmodels)
                Mu0 = Mu0';
            else
                error('Mu0 must be a column or row vector with the same number of dimensions as the number of states');
            end
            for s=1:nmodels
                 [X_p{s}(ind{s},1),W_p{s}(ind{s},ind{s},1)] = DecodingAlgorithms.PPDecode_predict(x0{s}(ind{s}), Pi0{s}(ind{s},ind{s}), Amat{s}(ind{s},ind{s},min(size(Amat{s},3),1)), Qmat{s}(:,:,min(size(Qmat{s},3))));
              
                if((estimateTarget==0 && ~isempty(yT{s})))               
                    invA= pinv(Amat{s}(:,:,min(size(Amat,3),1)));
                    ut{s}(:,1) = (Q{s}*pinv(PitT{s}(:,:,1)))*PhitT{s}(:,:,1)*(yT{s}-pinv(invA*PhitT{s}(:,:,1))*x0{s});
                    X_p{s}(ind{s},1) = X_p{s}(ind{s},1)+ut{s}(:,1);
                    W_p{s}(ind{s},ind{s},1) =W_p{s}(ind{s},ind{s},1) + (Q{s}*pinv(PitT{s}(:,:,1)))*A{s}*Pi0{s}*A{s}'*(Q{s}*pinv(PitT{s}(:,:,1)))';
                end
            end

    %            [~, S_est(1)] = max(MU_p(:,1)); %Most likely current state

                %State transition Probabilities must integrate to 1
                sumVal = sum(p_ij,2);
                if(any(sumVal~=1))
                    error('State Transition probability matrix must sum to 1 along each row');
                end
            %% 9 Steps
            % Filtering steps.
            for k = 1:(size(dN,2))

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 1 - p(s_k | H_k) = Sum(p(s_k|s_k-1)*p(s_k-1|H_k))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %MU(:,k)_p= [ p(s_k=1 | H_k)] is a vector of the probabilities
                %           [ p(s_k=2 | H_k)]
                %               .
                %               .
                %               .
                %           [ p(s_k=N | H_k)]
                % thus it is an prediction of the discrete state at time k given all
                % of the neural firing up to k-1 as summarized in H_k
                %
                % Whereas 
                % MU_u(:,k)=[ p(s_k=1 | H_k+1)] is a vector of the probabilities
                %           [ p(s_k=2 | H_k+1)]
                %               .
                %               .
                %               .
                %           [ p(s_k=N | H_k+1)]
                % The s suffix indicates that this is a "smoothed" estimate of
                % the state given the firing up to time k summarized in H_k+1
                if k==1
                    MU_p(:,k) = p_ij'*Mu0;          %state probability prediction equation
                else
                    MU_p(:,k) = p_ij'*MU_u(:,k-1);  
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 2 - p(s_k-1 | s_k, H_k)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % This is a matrix with i,j entry indicating the probability that 
                % s_k-1 = j given than s_k = i
                %
                % MU_p is the normalization factor. The first column of the
                % matrix of probabilities is:
                %
                % P(s_k-1=1 | s_k=1,H_k) ~ P(s_k=1|s_k-1=1,H_k)*P(s_k-1=1|H_k)
                % P(s_k-1=1 | s_k=2,H_k) ~ P(s_k=2|s_k-1=1,H_k)*P(s_k-1=1|H_k)
                % 
                % And the second columns ... etc
                %
                % 
                % P(s_k-1=2 | s_k=1,H_k) ~ P(s_k=1|s_k-1=2,H_k)*P(s_k-1=1|H_k)
                % P(s_k-1=2 | s_k=2,H_k) ~ P(s_k=2|s_k-1=2,H_k)*P(s_k-1=1|H_k)

                if(k==1)
                    p_ij_s= p_ij.*(Mu0*ones(1,nmodels));%.*(ones(nmodels,1)*(1./MU_p(:,k))');
                else
                    p_ij_s= p_ij.*(MU_u(:,k-1)*ones(1,nmodels));%.*(ones(nmodels,1)*(1./MU_p(:,k))');
                end
    %          
                 % To avoid any numerical issues with roundoff, we normalize to
                 % 1 again
                 normFact = repmat(sum(p_ij_s,1),[nmodels 1]); %Every column must sum to 1

                 p_ij_s = p_ij_s./normFact;
    %              for i=1:length(normFact)
    %                  if(normFact(i)~=0)
    %                     p_ij_s(:,i) = p_ij_s(:,i)./normFact(i);
    %                  else %reset all the states to be equally likely (each row must sum to 1)
    %                     p_ij_s(:,i) = 1/nmodels*ones(nmodels,1);
    %                  end
    %              end         
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Step 3 - Approximate p(x_k-1 | s_k, H_k) with Gaussian 
                % approximation to Mixtures of Gaussians
                % Calculate the mixed state mean for each filter
                % This will be the initial states for the update step of the
                % Point Process Adaptive Filter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    for j = 1:nmodels
                        for i = 1:nmodels
                            if(k>1)
                                X_s{j}(ind{i},k) = X_s{j}(ind{i},k) + X_u{i}(ind{i},k-1)*p_ij_s(i,j);
                            else
                                X_s{j}(ind{i},k) = X_s{j}(ind{i},k) + x0{i}(ind{i})*p_ij_s(i,j); 
                            end
                        end
                    end

                        % Calculate the mixed state covariance for each filter

                    for j = 1:nmodels
                        for i = 1:nmodels
                            if(k>1)
                                W_s{j}(ind{i},ind{i},k) = W_s{j}(ind{i},ind{i},k) + (W_u{i}(ind{i},ind{i},k-1) + (X_u{i}(ind{i},k-1)-X_s{j}(ind{i},k))*(X_u{i}(ind{i},k-1)-X_s{j}(ind{i},k))')*p_ij_s(i,j);
                            else
                                W_s{j}(ind{i},ind{i},k) = W_s{j}(ind{i},ind{i},k) + (Pi0{i}(ind{i},ind{i})+ (x0{i}(ind{i})-X_s{j}(ind{i},k))*(x0{i}(ind{i})-X_s{j}(ind{i},k))')*p_ij_s(i,j);
                            end
                        end
                    end

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 4 - Approximate p(x_k+1 |s_k+1,n_k+1,H_k+1)
               % Uses a bank of nmodel point process filters
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            k
               for s=1:nmodels

                   % Prediction Step
                   [X_p{s}(ind{s},k),W_p{s}(ind{s},ind{s},k)] = DecodingAlgorithms.PPDecode_predict(X_s{s}(ind{s},k), W_s{s}(ind{s},ind{s},k), Amat{s}(:,:,min(size(Amat,3),k)), Qmat{s}(:,:,min(size(Qmat{s},3))));

                   if(estimateTarget==0 && ~isempty(yT{s}))
                       if(k>1)
                            ut{s}(:,k) = (Q{s}*pinv(PitT{s}(:,:,k)))*PhitT{s}(:,:,k)*(yT{s}-pinv(PhitT{s}(:,:,k-1))*X_s{s}(ind{s},k));
                       else
                           invA = pinv(A{s}(:,:,min(size(A{s},3),1)));
                           ut{s}(:,k) = (Q{s}*pinv(PitT{s}(:,:,1)))*PhitT{s}(:,:,1)*(yT{s}-pinv(invA*PhitT{s}(:,:,1))*X_s{s}(ind{s},k));
                       end
                       X_p{s}(ind{s},k) = X_p{s}(ind{s},k)+ut{s}(:,k);
                       W_p{s}(ind{s},ind{s},k) =W_p{s}(ind{s},ind{s},k) + (Q{s}*pinv(PitT{s}(:,:,k)))*A{s}*W_s{s}(ind{s},ind{s},k)*A{s}'*(Q{s}*pinv(PitT{s}(:,:,k)))';

                    end

                   % Update Step
                   % Fold in the neural firing in the current time step
                   [X_u{s}(ind{s},k),W_u{s}(ind{s},ind{s},k),lambdaDeltaMat{s}(:,k)] = DecodingAlgorithms.PPDecode_update(X_p{s}(ind{s},k),squeeze(W_p{s}(ind{s},ind{s},k)),dN(:,1:k),lambdaCIFColl,binwidth,k);


               end

    %            
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 5 - p(n_k | s_k, H_k) using Laplace approximation
               % See General-purpose filter design for neural prosthetic devices.
               % Srinivasan L, Eden UT, Mitter SK, Brown EN.
               % J Neurophysiol. 2007 Oct;98(4):2456-75. Epub 2007 May 23.
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


               for s=1:nmodels

                 tempPdf = sqrt(det(W_u{s}(:,:,k)))./sqrt(det(W_p{s}(:,:,k)))*prod(exp(dN(:,k).*log(lambdaDeltaMat{s}(:,k))-lambdaDeltaMat{s}(:,k)));
                 pNGivenS(s,k) = tempPdf;
               end
               tempData = pNGivenS(:,k);
               tempData(isinf(tempData))=0;
               pNGivenS(:,k) = tempData;

               normFact = sum(pNGivenS(:,k));
               if(normFact~=0 && ~isnan(normFact))
                pNGivenS(:,k)=pNGivenS(:,k)./sum(pNGivenS(:,k));
               else

                   if(k>1)
                       pNGivenS(:,k) = pNGivenS(:,k-1);
                   else
                       pNGivenS(:,k) = .5*ones(nmodels,1);
                   end
               end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Step 6 - Calculate p(s_k | n_k, H_k) = p(s_k | H_k+1)
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               pSGivenN(:,k) = MU_p(:,k).*pNGivenS(:,k);



               %Normalization Factor
               normFact = sum(pSGivenN(:,k));
               if(normFact~=0 && ~isnan(normFact))
                pSGivenN(:,k) = pSGivenN(:,k)./sum(pSGivenN(:,k));
               else
                   if(k>1)
                        pSGivenN(:,k) = pSGivenN(:,k-1);
                   else
                        pSGivenN(:,k) = Mu0;
                   end

               end


               MU_u(:,k) = pSGivenN(:,k); %estimate of s_k given data up to k


               [~, S_est(k)] = max(MU_u(:,k)); %Most likely current state

               if(MinClassificationError==1)

                   s= S_est(k);
                   X(ind{s},k) = X_u{s}(ind{s},k);
                   W(ind{s},ind{s},k) = W_u{s}(ind{s},ind{s},k);

               else %Minimize the mean squared error

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   % Step 7 - Calculate p(x_k | n_k, H_k) - using gaussian
                   % approximation to mixture of gaussians 
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   for s=1:nmodels
                       X(ind{s},k) = X(ind{s},k)+MU_u(s,k)*X_u{s}(ind{s},k); 
                   end
                   for s=1:nmodels
                       W(ind{s},ind{s},k) =  W(ind{s},ind{s},k) +MU_u(s,k)*(W_u{s}(ind{s},ind{s},k) + (X_u{s}(ind{s},k)-X(ind{s},k))*(X_u{s}(ind{s},k)-X(ind{s},k))');
                   end
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               end
            end



         end   
        %% Kalman Filter 
        function [x_p, Pe_p, x_u, Pe_u,Gn,GnConvIter] = kalman_filter(A, C, Pv, Pw, Px0, x0,y, GnConv)
            %% DT Kalman Filter
            % This implements the DT Kalman filter for the system described by
            %
            % x(:,n+1) = A(:,:,n)x(:,n) + v(:,n)
            % y(:,n) = C(:,:,n)x(:,n) + w(:,n)
            %
            % where Pv(:,:,n), Pw(:,:,n) are the covariances of v(:,n) and w(:,n)
            % and Px0 is the initial state covariance.
            %
            % v(:,n), w(:,n), x(:,1) are assumed to be zero-mean.
            %
            % Return values are
            % x_p: state estimates given the past
            % Pe_p: error covariance estimates given the past
            % x_u: state updates given the data
            % Pe_u: error covariance updates given the data

            if(nargin<8||isempty(GnConv))
                GnConv = [];
            end
            N       = size(y,2); % number of time samples in the data
            x_p     = zeros( size(A,2), N+1 );
            x_u     = zeros( size(A,2), N );
            Pe_p    = zeros( size(A,2), size(A,2), N+1 );
            Gn      = zeros( size(A,2), size(C,1), N );
            Pe_u    = zeros( size(A,2), size(A,2), N );
            x_p(:,1)= x0;
            Pe_p(:,:,1) = Px0;

            for n=1:N
                [x_u(:,n),   Pe_u(:,:,n), Gn(:,:,n)]   = kalman_update( x_p(:,n), Pe_p(:,:,n), C(:,:,min(size(C,3),n)), Pw(:,:,min(size(Pw,3),n)), y(:,n),GnConv);
                [x_p(:,n+1), Pe_p(:,:,n+1)] = kalman_predict(x_u(:,n), Pe_u(:,:,n), A(:,:,min(size(A,3),n)), Pv(:,:,min(size(Pv,3),n)),GnConv);
                if(n>1 && isempty(GnConv))
                    diffGn = abs(Gn(:,:,n)-Gn(:,:,n-1));
                    mAbsdiffGn = max(max(diffGn));
                    if(mAbsdiffGn<1e-6)
                        GnConv=Gn(:,:,n);
                        GnConvIter = n;
                    else
                        GnConvIter=[];
                    end
                    
                end
            end

        


            %% Kalman Filter Update Equation
                function [x_u, Pe_u, G] = kalman_update(x_p, Pe_p, C, Pw, y, GnConv)
                % The Kalman update step that finds the state estimate based on new
                % data
                    if(nargin<6 || isempty(GnConv))
                        G       = (Pe_p * C')/(C * Pe_p * C' + Pw);
                    else 
                        G       = GnConv;
                    end
                    x_u     = x_p + G * (y - C * x_p);
                    Pe_u    = Pe_p - G * C * Pe_p;
                    Pe_u    = 0.5*(Pe_u + Pe_u');

%                     figure(10); subplot(1,3,1);imagesc(Pe_u); pause(0.005);
                end
            %% Kalman Filter Prediction Step
                function [x_p, Pe_p] = kalman_predict(x_u, Pe_u, A, Pv,GnConv)
                % The Kalman prediction step that implements the tracking system
                    x_p     = A * x_u;
                    if(isempty(GnConv))
                        Pe_p    = A * Pe_u * A' + Pv;
                    else
                        Pe_p    = Pe_u;
                    end
                    Pe_p    = 0.5*(Pe_p + Pe_p');
                end
        end
        %% Kalman Fixed-Interval Smoother
%         function  [x_pLag, Pe_pLag, x_uLag, Pe_uLag] = kalman_fixedIntervalSmoother(A, C, Pv, Pw, Px0, x0,y,lags)
%             %y should be zero mean gaussian
%             % only works for fixed A, C;
%             N       = size(y,2); % number of time samples in the data
%             nStates = size(A,2);
%             x_p  = zeros( (lags+1)*nStates, N+1 );
%             x_u  = zeros( (lags+1)*nStates, N );
%             Pe_p = zeros( (lags+1)*nStates, size(A,2), N+1 );
% %             GnLag   = zeros( (lags+1)*nStates, size(C,1), N );
%             Pe_u = zeros( (lags+1)*nStates, size(A,2), N );
%             x_p(1:nStates,1)= x0;
%             Pe_p(1:nStates,1:nStates,1) = Px0;
% 
%             for n=1:N
%                 [x_u(:,n),   Pe_u(:,:,n)]   = kalman_smootherUpdate( x_p(:,n), Pe_p(:,:,n), C, Pw, y(:,n),lags);
%                 Pe_pLag(:,:,n)= ((A')^-lags)*Pe_p((1:nStates)+(lags*nStates),(1:nStates),n); %W(n-N,n|n-1) -->W(n-N,n-N|n-1)
%                 %W(n-N,n-N|n-1) --> %W(n-N,n-N|n)
%                 Pe_uLag(:,:,n)= Pe_pLag(:,:,n)-Pe_p((1:nStates)+(lags*nStates),(1:nStates),n)*C'/(C*Pe_p((1:nStates),(1:nStates),n)*C'+Pw)*C*Pe_p((1:nStates)+(lags*nStates),(1:nStates),n)';             
%                 [x_p(:,n+1), Pe_p(:,:,n+1)] = kalman_smootherPredict(x_u(:,n), Pe_u(:,:,n), A, Pv,lags);
%                 %W(n-(N-1),n-(N-1)|n) <---- %W(n-N,n-N|n)
%                 Pe_pLag(:,:,n+1)=A*Pe_uLag(:,:,n)*A'+Pv;
%                 
%             end
%             offset = (lags*nStates);
%             x_pLag = x_p((1:nStates)+offset,:);
%             x_uLag = x_u((1:nStates)+offset,:);
%           
%              %% Kalman Filter Update Equation
%                 function [x_u, Pe_u] = kalman_smootherUpdate(x_p, Pe_p, C, Pw, y, lags)
%                 % The Kalman update step that finds the state estimate based on new
%                 % data
%                 nStates = size(C,2);
%                 I = eye(size(Pe_p(1:nStates,1:nStates)));
%                 tempMat=(C'/(C*Pe_p(1:nStates,1:nStates)*C'+Pw));
%                 tempVal=tempMat*(y-C*x_p(1:nStates));
%                 offset=0;
% 
%                 Pe_u    = Pe_p.*repmat((I- Pe_p((1:nStates)+offset,1:nStates)*tempMat*C*Pe_p((1:nStates),1:nStates)),lags+1,1);
%                 x_u     = x_p + Pe_p*tempVal;
% 
%                 Pe_u(1:nStates,1:nStates)    = 0.5*(Pe_u(1:nStates,1:nStates)  + Pe_u(1:nStates,1:nStates)');
%                     
%                    
%                    
%                 end
%             %% Kalman Filter Prediction Step
%                 function [x_p, Pe_p] = kalman_smootherPredict(x_u, Pe_u, A, Pv, lags)
%                 % The Kalman prediction step that implements the tracking system
%                 nStates = size(A,2);
% 
%                 Alag=zeros(length(x_u),length(x_u));
%                 Alag(1:nStates,1:nStates)=A;
%                 Alag((1:lags*nStates)+nStates,1:lags*nStates)=eye(lags*nStates,lags*nStates);
%                 Qlag=zeros(length(x_u),nStates);
%                 Qlag(1:nStates,1:nStates)=Pv;
%                 
%                 x_p = Alag*x_u;
%                 Pe_p =Alag*Pe_u*A' + Qlag;
%                 end
%                 
%         end
        %% Kalman Fixed-Interval Smoother
        function  [x_pLag, Pe_pLag, x_uLag, Pe_uLag] = kalman_fixedIntervalSmoother(A, C, Pv, Pw, Px0, x0,y,lags)
            %y should be zero mean gaussian
            N       = size(y,2);
            nStates = size(A,2);
            nObs = size(C,1);
            Alag = zeros((lags+1)*nStates,(lags+1)*nStates,N);
            Pvlag = zeros((lags+1)*nStates,(lags+1)*nStates,N);
            Clag = zeros(nObs,(lags+1)*nStates,N);
            Pwlag = zeros(nObs,nObs,N);
            x0lag = zeros(length(x0)*(lags+1),1);
            Px0lag = zeros((lags+1)*nStates,(lags+1)*nStates);
            Px0lag((1:nStates),(1:nStates))=Px0;
            x0lag(1:nStates,1)=x0;
            for n=1:N
                offset = 0;
                for i=1:(lags+1)
                    if(i==1)
                        Alag((1:nStates)+offset,(1:nStates)+offset,n)=A(:,:,min(size(A,3),n));
                        Pvlag((1:nStates)+offset,(1:nStates)+offset,n)=Pv(:,:,min(size(Pv,3),n));
                        Clag((1:nObs),(1:nStates)+offset,n)=C(:,:,min(size(C,3),n));
                        Pwlag((1:nObs),(1:nObs),n) = Pw(:,:,min(size(Pw,3),n));
                    else
                        Alag((1:nStates)+offset,(1:nStates)+(offset-nStates),n)=eye(nStates,nStates);
                        Pvlag((1:nStates)+offset,(1:nStates)+offset,n)=zeros(nStates,nStates);
                        Clag((1:nObs),(1:nStates)+offset,n)=zeros(nObs,nStates);
                    end
                    offset=offset+nStates;
                end
            end
           
            [x_p, Pe_p, x_u, Pe_u] = DecodingAlgorithms.kalman_filter(Alag, Clag, Pvlag, Pwlag, Px0lag, x0lag,y);

            x_pLag = x_p((lags*nStates+1):(lags+1)*nStates,:);
            Pe_pLag = Pe_p((lags*nStates+1):(lags+1)*nStates,(lags*nStates+1):(lags+1)*nStates,:);
            x_uLag = x_u((lags*nStates+1):(lags+1)*nStates,:);
            Pe_uLag = Pe_u((lags*nStates+1):(lags+1)*nStates,(lags*nStates+1):(lags+1)*nStates,:);
        end
        %% Kalman Smoother
        function [x_N, P_N,Ln] = kalman_smootherFromFiltered(A, x_p, Pe_p, x_u, Pe_u)
            N=size(x_u,2);

            x_N=zeros(size(x_u));
            P_N=zeros(size(Pe_u));
            Ln = zeros(size(P_N,1),size(P_N,2),size(P_N,3)-1);
            j=fliplr(1:N-1);
            x_N(:,N) = x_u(:,N);
            P_N(:,:,N) = Pe_u(:,:,N);
%             LnConv = [];
            for n=j
%                 if(n<round(N/100) || N<10000)
                    Ln(:,:,n)=Pe_u(:,:,n)*A(:,:,min(size(A,3),n))'/Pe_p(:,:,n+1);
%                 elseif(~isempty(LnConv))
%                     Ln(:,:,n)=LnConv;
%                 else
%                     Ln(:,:,n)=Pe_u(:,:,n)*A(:,:,min(size(A,3),n))'/Pe_p(:,:,n+1);
%                 end
                x_N(:,n) = x_u(:,n)+Ln(:,:,n)*(x_N(:,n+1)-x_p(:,n+1));
                P_N(:,:,n)=Pe_u(:,:,n)+Ln(:,:,n)*(P_N(:,:,n+1)-Pe_p(:,:,n+1))*Ln(:,:,n)';
                P_N(:,:,n) = 0.5*(P_N(:,:,n)+P_N(:,:,n)');
%                 if(n<(N-1) && isempty(LnConv))
%                     diffLn = abs(Ln(:,:,n)-Ln(:,:,n+1));
%                     mAbsdiffLn = max(max(diffLn));
%                     if(mAbsdiffLn<1e-6)
%                         LnConv=Ln(:,:,n);
%                         LnConvIter = n;
%                     end
%                 end
            end    
   
        
         end
        function [x_N, P_N,Ln,x_p, Pe_p, x_u, Pe_u] = kalman_smoother(A, C, Pv, Pw, Px0, x0, y)
            %% kalman smoother
            N=size(y,2);
            [x_p, Pe_p, x_u, Pe_u] = kalman_filter(A, C, Pv, Pw, Px0, x0, y);

            x_N=zeros(size(x_u));
            P_N=zeros(size(Pe_u));
            Ln = zeros(size(P_N,1),size(P_N,2),size(P_N,3)-1);
            j=fliplr(1:N-1);
            x_N(:,N) = x_u(:,N);
            P_N(:,:,N) = Pe_u(:,:,N);
%             LnConv = [];
            for n=j
%                 if(n<round(N/100)|| N<10000)
                    Ln(:,:,n)=Pe_u(:,:,n)*A(:,:,min(size(A,3),n))'/Pe_p(:,:,n+1);
%                 elseif(~isempty(LnConv))
%                     Ln(:,:,n)=LnConv;
%                 else
%                     Ln(:,:,n)=Pe_u(:,:,n)*A(:,:,min(size(A,3),n))'/Pe_p(:,:,n+1);
%                 end
                x_N(:,n) = x_u(:,n)+Ln(:,:,n)*(x_N(:,n+1)-x_p(:,n+1));
                P_N(:,:,n)=Pe_u(:,:,n)+Ln(:,:,n)*(P_N(:,:,n+1)-Pe_p(:,:,n+1))*Ln(:,:,n)';
                P_N(:,:,n) = 0.5*(P_N(:,:,n)+P_N(:,:,n)');
%                 if(n<(N-1) && isempty(LnConv))
%                     diffLn = abs(Ln(:,:,n)-Ln(:,:,n+1));
%                     mAbsdiffLn = max(max(diffLn));
%                     if(mAbsdiffLn<1e-6)
%                         LnConv=Ln(:,:,n);
%                         LnConvIter = n;
%                     end
%                 end
            end
        end

        %% Functions for Point Process State Space Expectation Maximization 
        % PPSS_EMFB implements a forward and backward PPSS_EM. Because the way that the algorithm is setup,
        % we can analyze the data from the first trial to the last (forward) and from
        % the last trial to the first (backward). This approach yields
        % better estimates of the underlying firing rates
        function [xKFinal,WKFinal, WkuFinal,Qhat,gammahat,fitResults,stimulus,stimCIs,logll,QhatAll,gammahatAll,nIter]=PPSS_EMFB(A,Q0,x0,dN,fitType,delta,gamma0,windowTimes, numBasis,neuronName)
    %         if(nargin<10 || isempty(neuronName))
    %             neuronName = 1;
    %         end
            dLikelihood(1)=inf;
            if(numel(Q0)==length(Q0)^2)
                Q0=diag(Q0); %turn Q into a vector
            end

            Qhat=Q0;
            gammahat=gamma0;
            xK0=x0;
            cnt=1; tol=1e-2; maxIter=2e3;
            tolAbs = 1e-3;
            tolRel = 1e-3;
            llTol  = 1e-3;
            stoppingCriteria=0;

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    HkAll{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    HkAll{k} = 0;
                end
                gamma0=0;
            end


            HkAllR=HkAll(end:-1:1);
    %         if(~isempty(windowTimes))
    %             histObj = History(windowTimes,minTime,maxTime);
    %             for k=K:-11:1
    %                 nstr{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    %                 nstr{k}.setMinTime(minTime);
    %                 nstr{k}.setMaxTime(maxTime);
    %                 HkAllR{k} = histObj.computeHistory(nstr{k}).dataToMatrix;
    %             end
    %         else
    %             for k=1:K
    %                 HkAllR{k} = 0;
    %             end
    %             gammahat=0;
    %         end

            while(stoppingCriteria~=1 && cnt<maxIter)
                display('EMFB: Forward EM');
                [xK,WK, Wku,Qhat(:,cnt+1),gammahat(cnt+1,:),logll(cnt),~,~,nIter1,negLL]=DecodingAlgorithms.PPSS_EM(A,Qhat(:,cnt),xK0,dN,fitType,delta,gammahat(cnt,:),windowTimes, numBasis,HkAll);    
                if(~negLL)
                    display('EMFB: Backward EM');
                    [xKR,~, ~,QhatR(:,cnt+1),gammahatR(cnt+1,:),logllR(cnt),~,~,nIter2,negLL]=DecodingAlgorithms.PPSS_EM(A,Qhat(:,cnt+1),xK(:,end),flipud(dN),fitType,delta,gammahat(cnt+1,:),windowTimes, numBasis,HkAllR);
                    if(~negLL)
                        display('EMFB: Forward EM');

                        [xK2,WK2, Wku2,Qhat2,gammahat2,logll2,~,~,nIter3,negLL2]=DecodingAlgorithms.PPSS_EM(A,QhatR(:,cnt+1),xKR(:,end),dN,fitType,delta,gammahatR(cnt+1,:),windowTimes, numBasis,HkAll);
                        if(~negLL2)
                            xK=xK2;
                            WK=WK2;
                            Wku=Wku2;
                            Qhat(:,cnt+1) = Qhat2;
                            gammahat(cnt+1,:) = gammahat2;
                            logll(cnt) = logll2;

                        end
                    end
                end


                xK0=xK(:,1);
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(logll(cnt)-logll(cnt-1));%./abs(logll(cnt-1));
                end
                cnt=cnt+1;

    %             figure(1)
    %         
    %             subplot(1,2,1); surf(xK);
    %             subplot(1,2,2); plot(logll); ylabel('Log Likelihood');

                dQvals = abs(sqrt(Qhat(:,cnt))-sqrt(Qhat(:,cnt-1)));
                dGamma = abs(gammahat(cnt,:)-gammahat(cnt-1,:));
                dMax = max([dQvals',dGamma]);

                dQRel = max(abs(dQvals./sqrt(Qhat(:,cnt-1))));
                dGammaRel = max(abs(dGamma./gammahat(cnt-1,:)));
                dMaxRel = max([dQRel,dGammaRel]);
    %             dMax
    %             dMaxRel
                if(dMax<tolAbs && dMaxRel<tolRel)
                    stoppingCriteria=1;
                    display(['EMFB converged at iteration:' num2str(cnt) ' b/c change in params was within criteria']);
                end
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
                    stoppingCriteria=1;
                    display(['EMFB stopped at iteration:' num2str(cnt) ' b/c change in likelihood was negative']);
                end

            end

            maxLLIndex = find(logll == max(logll),1,'first');
            if(maxLLIndex==1)
                maxLLIndex=cnt-1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
            end

            xKFinal = xK;
            x0Final=xK(:,1);
            WKFinal = WK;
            WkuFinal = Wku;
            QhatAll =Qhat(:,1:maxLLIndex+1);
            Qhat = Qhat(:,maxLLIndex+1);
            gammahatAll =gammahat(1:maxLLIndex+1);
            gammahat = gammahat(maxLLIndex+1,:);
            logll = logll(maxLLIndex);

            K=size(dN,1);
            SumXkTermsFinal = diag(Qhat(:,:,end))*K;
            logllFinal=logll(end);
            McInfo=100;
            McCI = 3000;

            nIter = [];%[nIter1,nIter2,nIter3];
  
            
            K   = size(dN,1); 
            R=size(xK,1);
            logllobs = logll+R*K*log(2*pi)+K/2*log(det(diag(Qhat)))+ 1/2*trace(diag(Qhat)\SumXkTermsFinal);

            InfoMat = DecodingAlgorithms.estimateInfoMat(fitType,dN,HkAll,A,x0Final,xKFinal,WKFinal,WkuFinal,Qhat,gammahat,windowTimes,SumXkTermsFinal,delta,McInfo);
            fitResults = DecodingAlgorithms.prepareEMResults(fitType,neuronName,dN,HkAll,xKFinal,WKFinal,Qhat,gammahat,windowTimes,delta,InfoMat,logllobs);
            [stimCIs, stimulus] = DecodingAlgorithms.ComputeStimulusCIs(fitType,xKFinal,WkuFinal,delta,McCI);
    %             

        end
        function [xKFinal,WKFinal, WkuFinal,Qhat,gammahat,logll,QhatAll,gammahatAll,nIter,negLL]=PPSS_EM(A,Q0,x0,dN,fitType,delta,gamma0,windowTimes, numBasis,Hk)
            if(nargin<9 || isempty(numBasis))
                numBasis = 20;
            end
            if(nargin<8 || isempty(windowTimes))
                if(isempty(gamma0))
                    windowTimes =[];
                else
    %                 numWindows =length(gamma0)+1; 
                    windowTimes = 0:delta:(length(gamma0)+1)*delta;
                end
            end
            if(nargin<7)
                gamma0=[];
            end
            if(nargin<6 || isempty(delta))
                delta = .001;
            end
            if(nargin<5)
                fitType = 'poisson';
            end


            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            K=size(dN,1);




    %         tol = 1e-3; %absolute change;
            tolAbs = 1e-3;
            tolRel = 1e-3;
            llTol  = 1e-3;
            cnt=1;

            maxIter = 100;

            if(numel(Q0)==length(Q0)^2)
                Q0=diag(Q0); %turn Q into a vector
            end
               numToKeep=10;
            Qhat = zeros(length(Q0),numToKeep);
            Qhat(:,1)=Q0;
            gammahat=zeros(numToKeep,length(gamma0));
            gammahat(1,:)=gamma0;
%             QhatNew=Q0;
%             gammahatNew(1,:)=gamma0;
            cnt=1;
            dLikelihood(1)=inf;
    %         logll(1)=-inf;
            x0hat = x0;
            negLL=0;
         
            %Forward EM
            stoppingCriteria =0;
%             logllNew= -inf;
            while(stoppingCriteria~=1 && cnt<=maxIter)
                 storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
                 storeIndP1= mod(cnt,numToKeep)+1;
                 storeIndM1= mod(cnt-2,numToKeep)+1;
                [xK{storeInd},WK{storeInd},Wku{storeInd},logll(cnt),SumXkTerms,sumPPll]= ...
                    DecodingAlgorithms.PPSS_EStep(A,Qhat(:,storeInd),x0hat,dN,Hk,fitType,delta,gammahat(storeInd,:),numBasis);
             
                [Qhat(:,storeIndP1),gammahat(storeIndP1,:)] = DecodingAlgorithms.PPSS_MStep(dN,Hk,fitType,xK{storeInd},WK{storeInd},gammahat(storeInd,:),delta,SumXkTerms,windowTimes);
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(logll(cnt)-logll(cnt-1));%./abs(logll(cnt-1));
                end

                if(mod(cnt,25)==0)
                    figure(1);
                    subplot(1,2,1); surf(xK{storeInd});
                    subplot(1,2,2); plot(logll); ylabel('Log Likelihood');
                end
                
                dQvals = abs(sqrt(Qhat(:,storeInd))-sqrt(Qhat(:,storeIndM1)));
                dGamma = abs(gammahat(storeInd,:)-gammahat(storeIndM1,:));
                dMax = max([dQvals',dGamma]);

                dQRel = max(abs(dQvals./sqrt(Qhat(:,storeIndM1))));
                dGammaRel = max(abs(dGamma./gammahat(storeIndM1,:)));
                dMaxRel = max([dQRel,dGammaRel]);

             
                cnt=(cnt+1);
                if(dMax<tolAbs && dMaxRel<tolRel)
                    stoppingCriteria=1;
                    display(['         EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
                    negLL=0;
                end
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
                    stoppingCriteria=1;
                    display(['         EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
                    negLL=1;
                end
                

            end


            maxLLIndex  = find(logll == max(logll),1,'first');
            maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
            if(maxLLIndex==1)
%                 maxLLIndex=cnt-1;
                maxLLIndex =1;
                maxLLIndMod = 1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
               maxLLIndMod = 1;
%             else
%                maxLLIndMod = mod(maxLLIndex,numToKeep); 
               
            end
            nIter   = cnt-1;  
%             maxLLIndMod
            xKFinal = xK{maxLLIndMod};
            WKFinal = WK{maxLLIndMod};
            WkuFinal = Wku{maxLLIndMod};
            QhatAll =Qhat(:,1:maxLLIndMod);
            Qhat = Qhat(:,maxLLIndMod);
            gammahatAll =gammahat(1:maxLLIndMod);
            gammahat = gammahat(maxLLIndMod,:);
            logll = logll(maxLLIndex);
           
        end
        
        % Subroutines for the PPSS_EM algorithm
        function [x_K,W_K,Wku,logll,sumXkTerms,sumPPll]=PPSS_EStep(A,Q,x0,dN,HkAll,fitType,delta,gamma,numBasis)


             minTime=0;
             maxTime=(size(dN,2)-1)*delta;


             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end
            if(numel(Q)==length(Q))
                Q=diag(Q); %turn Q into a diagonal matrix
            end
            [K,N]   = size(dN); 
            R=size(basisMat,2);

            x_p     = zeros( size(A,2), K );
            x_u     = zeros( size(A,2), K );
            W_p    = zeros( size(A,2),size(A,2), K);
            W_u    = zeros( size(A,2),size(A,2), K );



            for k=1:K

                if(k==1)
                    x_p(:,k)     = A * x0;
                    W_p(:,:,k)   = Q;
                else
                    x_p(:,k)     = A * x_u(:,k-1);
                    W_p(:,:,k)   = A * W_u(:,:,k-1) * A' + Q;
                end

                 sumValVec=zeros(size(W_p,1),1);
                 sumValMat=zeros(size(W_p,2),size(W_p,2));



                if(strcmp(fitType,'poisson'))
                    Hk=HkAll{k};
                    Wk = basisMat*diag(W_p(:,:,k));
                    stimK=basisMat*x_p(:,k);

                    histEffect=exp(gamma*Hk')';
                    stimEffect=exp(stimK);
                    lambdaDelta =stimEffect.*histEffect;
                    GradLogLD =basisMat;
                    JacobianLogLD = zeros(R,R);
                    GradLD = basisMat.*repmat(lambdaDelta,[1 R]);

                    sumValVec = GradLogLD'*dN(k,:)' - diag(GradLD'*basisMat);
                    sumValMat = GradLD'*basisMat;


                elseif(strcmp(fitType,'binomial'))
                    Hk=HkAll{k};
                    Wk = basisMat*diag(W_p(:,:,k));
                    stimK=basisMat*x_p(:,k);

                    lambdaDelta=exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));  
                    GradLogLD =basisMat.*(repmat(1-lambdaDelta,[1 R]));
                    JacobianLogLD = basisMat.*repmat(lambdaDelta.*(-1+lambdaDelta),[1 R]);
                    GradLD = basisMat.*(repmat(lambdaDelta.*(1-lambdaDelta),[1 R]));
                    JacobianLD = basisMat.*(repmat(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta.^2),[1 R]));

                    sumValVec = GradLogLD'*dN(k,:)' - diag(GradLD'*basisMat);
                    sumValMat = -diag(JacobianLogLD'*dN(k,:)')+ JacobianLD'*basisMat;

                end  


%                  invW_u             = pinv(W_p(:,:,k))+ sumValMat;
%                  W_u(:,:,k)       = pinv(invW_u);% +100*diag(eps*rand(size(W_p,1),1));
%                  
                 

                 invW_u             = eye(size(W_p(:,:,k)))/W_p(:,:,k)+ sumValMat;
                 W_u(:,:,k)       = eye(size(invW_u))/invW_u;% +100*diag(eps*rand(size(W_p,1),1));


                 % Maintain Positive Definiteness
                % Make sure eigenvalues are positive
                [vec,val]=eig(W_u(:,:,k) ); val(val<=0)=eps;
                W_u(:,:,k) =vec*val*vec';
                x_u(:,k)      = x_p(:,k)  + W_u(:,:,k)*(sumValVec);

            end

            [x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u);


            Wku=zeros(R,R,K,K);
            Tk = zeros(R,R,K-1);
            for k=1:K
                Wku(:,:,k,k)=W_K(:,:,k);
            end

            for u=K:-1:2
                for k=(u-1):-1:1
                    Tk(:,:,k)=A;
%                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
                    Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
                    Wku(:,:,u,k)=Wku(:,:,k,u);
                end
            end
   

            %All terms
            Sxkxkp1 = zeros(R,R);
            Sxkp1xkp1 = zeros(R,R);
            Sxkxk = zeros(R,R);
            for k=1:K-1
%                Sxkxkp1 = Sxkxkp1+x_u(:,k)*x_K(:,k+1)'+ ...
%                    Lk(:,:,k)*(W_K(:,:,k+1)+(x_K(:,k+1)-x_p(:,k+1))*x_K(:,k+1)');
               Sxkxkp1 = Sxkxkp1+Wku(:,:,k,k+1)+x_K(:,k)*x_K(:,k+1)';
               Sxkp1xkp1 = Sxkp1xkp1+W_K(:,:,k+1)+x_K(:,k+1)*x_K(:,k+1)';
               Sxkxk = Sxkxk+W_K(:,:,k)+x_K(:,k)*x_K(:,k)';

            end

            sumXkTerms =  Sxkp1xkp1-A*Sxkxkp1-Sxkxkp1'*A'+A*Sxkxk*A'+ ...
                          W_K(:,:,1)+x_K(:,1)*x_K(:,1)' + ... %expected value of xK(1)^2
                          -A*x0*x_K(:,1)' -x_K(:,1)*x0'*A' +A*(x0*x0')*A';


            if(strcmp(fitType,'poisson'))
                sumPPll=0;
                for k=1:K
                    Hk=HkAll{k};
                    Wk = basisMat*diag(W_K(:,:,k));
                    stimK=basisMat*x_K(:,k);
                    histEffect=exp(gamma*Hk')';
                    stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
        %             stimEffect=exp(stimK  + Wk*0.5);
                    ExplambdaDelta =stimEffect.*histEffect;
                    ExplogLD = (stimK + (gamma*Hk')');
                    sumPPll=sum(dN(k,:)'.*ExplogLD - ExplambdaDelta);

                end
            elseif(strcmp(fitType,'binomial'))

                sumPPll=0;
                for k=1:K
                    Hk=HkAll{k};
                    Wk = basisMat*diag(W_K(:,:,k));
                    stimK=basisMat*x_K(:,k);
                    lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                    ExplambdaDelta=lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;  
            %       logLD = stimK+(gamma*Hk')' - log(1+lambdaDelta)
                    ExplogLD = stimK+(gamma*Hk')' - log(1+exp(stimK+(gamma*Hk')')) -Wk.*(lambdaDelta).*(1-lambdaDelta)*.5;
                    %E(f(x)]=f(x_hat) + 1/2sigma_x^2 * d^2/dx*f(x_hat)
                    %This is applied to log(1+exp(x_K))
            %         
                    sumPPll=sum(dN(k,:)'.*ExplogLD - ExplambdaDelta);
                end

            end
            R=numBasis;

            logll = -R*K*log(2*pi)-K/2*log(det(Q))  + sumPPll - 1/2*trace(pinv(Q)*sumXkTerms);


        end
        function [Qhat,gamma_new] = PPSS_MStep(dN,HkAll,fitType,x_K,W_K,gamma, delta,sumXkTerms,windowTimes)
             K=size(dN,1);
             N=size(dN,2);


             sumQ =  diag(diag(sumXkTerms));
             Qhat = sumQ*(1/K);

             [vec,val]=eig(Qhat); val(val<=0)=0.00000001;
             Qhat =vec*val*vec';
             Qhat = (diag(Qhat));


             minTime=0;
             maxTime=(size(dN,2)-1)*delta;

             numBasis = size(x_K,1);
             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end




            gamma_new = gamma;

            if(~isempty(windowTimes) && all(gamma_new~=0))
                converged=0;
                iter = 1;
                maxIter=300;
                while(~converged && iter<maxIter)
        %         disp(['      - Newton-Raphson alg. iter #',num2str(iter)])
                    if(strcmp(fitType,'poisson'))

                            gradQ=zeros(size(gamma_new,2),1);
                            jacQ =zeros(size(gamma_new,2),size(gamma_new,2));
                        for k=1:K
                            Hk=HkAll{k};

                            Wk = basisMat*diag(W_K(:,:,k));
                            stimK=basisMat*(x_K(:,k));
                            histEffect=exp(gamma_new*Hk')';
                            stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
    %                         stimEffect=exp(stimK+Wk*0.5);
                            lambdaDelta = stimEffect.*histEffect;

                            gradQ = gradQ + Hk'*dN(k,:)' - Hk'*lambdaDelta;
    %                         jacQ  = jacQ  - diag(diag((Hk.*repmat(lambdaDelta,[1 size(Hk,2)]))'*Hk));
                            jacQ  = jacQ  - (Hk.*repmat(lambdaDelta,[1 size(Hk,2)]))'*Hk;
                        end



                    elseif(strcmp(fitType,'binomial'))
                            gradQ=zeros(size(gamma_new,2),1);
                            jacQ =zeros(size(gamma_new,2),size(gamma_new,2));
                         for k=1:K
                            Hk=HkAll{k};

                            Wk = basisMat*diag(W_K(:,:,k));
                            stimK=basisMat*(x_K(:,k));

                            histEffect=exp(gamma_new*Hk')';
                            stimEffect=exp(stimK);
    %                         stimEffect=exp(stimK+Wk*0.5);
                            C = stimEffect.*histEffect;
                            M = 1./C;
                            lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                            ExpLambdaDelta = lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;
                            ExpLDSquaredTimesInvExp = (lambdaDelta).^2.*1./C;
                            ExpLDCubedTimesInvExpSquared = (lambdaDelta).^3.*M.^2 +Wk/2.*(3.*M.^4.*lambdaDelta.^3+12.*lambdaDelta.^3.*M.^3-12.*M.^4.*lambdaDelta.^4);

    %                         ExpLambdaDeltaTimesExp = C.*lambdaDelta + (2.*C.*lambdaDelta-3*C.*lambdaDelta.*lambdaDelta).*Wk/2;
    %                         ExpLambdaDeltaTimesExpSquared = C.^2.*lambdaDelta + (7.*C.^2.*lambdaDelta-5*C.^2.*lambdaDelta.*lambdaDelta).*Wk/2;

    %                         lambdaDelta = C./(1+C);

    %                         gradQ = gradQ + (Hk.*repmat(1-lambdaDelta,[1,size(Hk,2)]))'*dN(k,:)' ...
    %                                       - (Hk.*repmat(C,[1,size(Hk,2)]))'*lambdaDelta;
    %                         jacQ  = jacQ  - (Hk.*repmat(C.*lambdaDelta.*dN(k,:)',[1,size(Hk,2)]))'*Hk ...
    %                                       - (Hk.*repmat(lambdaDelta,[1,size(Hk,2)]))'*Hk ...
    %                                       - (Hk.*repmat(C.^2.*lambdaDelta,[1,size(Hk,2)]))'*Hk;
                            gradQ = gradQ + (Hk.*repmat(1-ExpLambdaDelta,[1,size(Hk,2)]))'*dN(k,:)' ...
                                          - (Hk.*repmat(ExpLDSquaredTimesInvExp./lambdaDelta,[1,size(Hk,2)]))'*lambdaDelta;
                            jacQ  = jacQ  - (Hk.*repmat(ExpLDSquaredTimesInvExp.*dN(k,:)',[1,size(Hk,2)]))'*Hk ...
                                          - (Hk.*repmat(ExpLDSquaredTimesInvExp,[1,size(Hk,2)]))'*Hk ...
                                          - (Hk.*repmat(2*ExpLDCubedTimesInvExpSquared,[1,size(Hk,2)]))'*Hk;

                         end


                    end

                    gamma_newTemp = (gamma_new'-pinv(jacQ)*gradQ)';
                    if(any(isnan(gamma_newTemp)))
                        gamma_newTemp = gamma_new;
    %                     gradQ=max(gradQ,-10);
    %                     gradQ=min(gradQ,10);
    %                     gamma_newTemp = (gamma_new' - jacQ\gradQ)';
    %                     if(any(isnan(gamma_newTemp)))
    %                         if(isinf(gamma_new))
    %                             gamma_newTemp(isinf(gamma_new))=-5;
    %                         else
    %                             gamma_newTemp=gamma_new;    
    %                         end
    %                         
    %                     end
    %                 elseif(abs(gamma_newTemp)>1e1)
    %                     gamma_newTemp = sign(gamma_newTemp)*1e1;
                    end
                    mabsDiff = max(abs(gamma_newTemp - gamma_new));
                    if(mabsDiff<10^-2)
                        converged=1;
                    end
                    gamma_new=gamma_newTemp;
                    iter=iter+1;
                end
                %Keep gamma from getting too large since this effect is
                %exponentiated
                gamma_new(gamma_new>1e2)=1e1;
                gamma_new(gamma_new<-1e2)=-1e1;
            end

    %          pause;
        end
        function fitResults=prepareEMResults(fitType,neuronNumber,dN,HkAll,xK,WK,Q,gamma,windowTimes,delta,informationMatrix,logll)


            [numBasis, K] =size(xK);
            SE = sqrt(abs(diag(inv(informationMatrix))));
            xKbeta = reshape(xK,[numel(xK) 1]);
            seXK=[];
            for k=1:K
                seXK   = [seXK; sqrt(diag(WK(:,:,k)))];
            end
            statsStruct.beta=[xKbeta;(Q(:,end));gamma(end,:)'];
            statsStruct.se  =[seXK;SE];
            covarianceLabels = cell(1,numBasis);
            for r=1:numBasis
                if(r<10)
                    covarianceLabels{r} =  ['Q0' num2str(r)];
                else
                    covarianceLabels{r} =  ['Q' num2str(r)];
                end
            end

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end

            nst = cell(1,K);
            if(~isempty(windowTimes))
                histObj{1} = History(windowTimes,minTime,maxTime);
            else
                histObj{1} = [];
            end

            if(isnumeric(neuronNumber))
                name=num2str(neuronNumber);
                if(neuronNumber>0 && neuronNumber<10)
                    name = strcat(num2str(0),name);
                end
                name = ['N' name];  
            else
                name = neuronNumber;
            end

            for k=1:K
                nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta,name);
                nst{k}.setMinTime(minTime);
                nst{k}.setMaxTime(maxTime);

            end

            nCopy = nstColl(nst);
            nCopy = nCopy.toSpikeTrain;
            lambdaData=[];
            cnt=1;

            for k=1:K
                Hk=HkAll{k};
                stimK=basisMat*xK(:,k);


                if(strcmp(fitType,'poisson'))
                    histEffect=exp(gamma(end,:)*Hk')';
                    stimEffect=exp(stimK);
                    lambdaDelta = histEffect.*stimEffect;
                    lambdaData = [lambdaData;lambdaDelta/delta];
                elseif(strcmp(fitType,'binomial'))
                    histEffect=exp(gamma(end,:)*Hk')';
                    stimEffect=exp(stimK);
                    lambdaDelta = histEffect.*stimEffect;
                    lambdaDelta = lambdaDelta./(1+lambdaDelta);
                    lambdaData = [lambdaData;lambdaDelta/delta];
                end


                for r=1:numBasis
                        if(r<10)
                            otherLabels{cnt} = ['b0' num2str(r) '_{' num2str(k) '}']; 
                        else
                            otherLabels{cnt} = ['b' num2str(r) '_{' num2str(k) '}'];
                        end
                        cnt=cnt+1;
                end
            end

            lambdaTime = minTime:delta:(length(lambdaData)-1)*delta;
            nCopy.setMaxTime(max(lambdaTime));
            nCopy.setMinTime(min(lambdaTime));

            numLabels = length(otherLabels);
            if(~isempty(windowTimes))
                histLabels  = histObj{1}.computeHistory(nst{1}).getCovLabelsFromMask;
            else
                histLabels = [];
            end
            otherLabels((numLabels+1):(numLabels+length(covarianceLabels)))=covarianceLabels;
            numLabels = length(otherLabels);

            tc{1} = TrialConfig(otherLabels,sampleRate,histObj,[]); 
            numBasisStr=num2str(numBasis);
            numHistStr = num2str(length(windowTimes)-1);
            if(~isempty(histObj))
                tc{1}.setName(['SSGLM(N_{b}=', numBasisStr,')+Hist(N_{h}=' ,numHistStr,')']);
            else
                tc{1}.setName(['SSGLM(N_{b}=', numBasisStr,')']);
            end
            configColl= ConfigColl(tc);


            otherLabels((numLabels+1):(numLabels+length(histLabels)))=histLabels;




            labels{1}  = otherLabels; % Labels change depending on presence/absense of History or ensCovHist
            if(~isempty(windowTimes))
                numHist{1} = length(histObj{1}.windowTimes)-1;
            else 
                numHist{1}=[];
            end

            ensHistObj{1} = [];
            lambdaIndexStr=1;
            lambda=Covariate(lambdaTime,lambdaData,...
                           '\Lambda(t)','time',...
                           's','Hz',strcat('\lambda_{',lambdaIndexStr,'}'));


            AIC = 2*length(otherLabels)-2*logll;
            BIC = -2*logll+length(otherLabels)*log(length(lambdaData));

            dev=-2*logll;
            b{1} = statsStruct.beta;
            stats{1} = statsStruct;

            distrib{1} =fitType;
            currSpikes=nst;%nspikeColl.getNST(tObj.getNeuronIndFromName(neuronNames));
            for n=1:length(currSpikes)
                currSpikes{n} = currSpikes{n}.nstCopy;
                currSpikes{n}.setName(nCopy.name);
            end
            XvalData{1} = [];
            XvalTime{1} = [];
            spikeTraining = currSpikes;


            fitResults=FitResult(spikeTraining,labels,numHist,histObj,ensHistObj,lambda,b, dev, stats,AIC,BIC,configColl,XvalData,XvalTime,distrib);
            DTCorrection=1;
            makePlot=0;
            Analysis.KSPlot(fitResults,DTCorrection,makePlot);
            Analysis.plotInvGausTrans(fitResults,makePlot);
            Analysis.plotFitResidual(fitResults,[],makePlot); 
        end
        function [CIs, stimulus]  = ComputeStimulusCIs(fitType,xK,Wku,delta,Mc,alphaVal)
            if(nargin<6 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<5 ||isempty(Mc))
                Mc=3000;
            end
            [numBasis,K]=size(xK);


           for r=1:numBasis  
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
    %                 stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))/delta;
                    if(strcmp(fitType,'poisson'))
                        stimulusDraw(r,:,c) =  exp(xKDraw(r,:,c))/delta;
                    elseif(strcmp(fitType,'binomial'))
                        stimulusDraw(r,:,c) = exp(xKDraw(r,:,c))./(1+exp(xKDraw(r,:,c)))/delta;
                    end
                end
           end

           CIs = zeros(size(xK,1),size(xK,2),2);
           for r=1:numBasis
               for k=1:K
                   [f,x] = ecdf(squeeze(stimulusDraw(r,k,:)));
                    CIs(r,k,1) = x(find(f<alphaVal/2,1,'last'));
                    CIs(r,k,2) = x(find(f>(1-(alphaVal/2)),1,'first'));
               end
           end

           if(nargout==2)
               if(strcmp(fitType,'poisson'))
                    stimulus =  exp(xK)/delta;
               elseif(strcmp(fitType,'binomial'))
                    stimulus = exp(xK)./(1+exp(xK))/delta;
               end
           end


        end
        function InfoMatrix=estimateInfoMat(fitType,dN,HkAll,A,x0,xK,WK,Wku,Q,gamma,windowTimes,SumXkTerms,delta,Mc)
            if(nargin<14)
                Mc=500;
            end

            [K,N]=size(dN);
            if(~isempty(windowTimes))
                J=max(size(gamma(end,:)));
            else
                J=0;
            end

            R=size(Q,1);
            numBasis = R;

            % The complete data information matrix
            Ic=zeros(J+R,J+R);
            Q=(diag(Q)); % Make sure Q is diagonal matrix


            X=((SumXkTerms));
            Ic(1:R,1:R) = K/2*eye(size(Q))/Q^2 +X'/Q^3;


            % Compute information of history terms
            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
    %         nst = cell(1,K);
    %         if(~isempty(windowTimes))
    %             histObj = History(windowTimes,minTime,maxTime);
    %             for k=1:K
    %                 nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
    %                 nst{k}.setMinTime(minTime);
    %                 nst{k}.setMaxTime(maxTime);
    %                 Hn{k} = histObj.computeHistory(nst{k}).dataToMatrix;
    %             end
    %         else
    %             for k=1:K
    %                 Hn{k} = 0;
    %             end
    %             gamma=0;
    %         end

             if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
             end

            jacQ =zeros(size(gamma,2),size(gamma,2));
            if(strcmp(fitType,'poisson'))
                for k=1:K
                    Hk=HkAll{k};

                    Wk = basisMat*diag(WK(:,:,k));
                    stimK=basisMat*(xK(:,k));
                    histEffect=exp(gamma*Hk')';
                    stimEffect=exp(stimK)+exp(stimK)/2.*Wk;
                    lambdaDelta = stimEffect.*histEffect;

                    jacQ  = jacQ  - (Hk.*repmat(lambdaDelta,[1 size(Hk,2)]))'*Hk;
                end

             elseif(strcmp(fitType,'binomial'))
                 for k=1:K
                    Hk=HkAll{k};
                    Wk = basisMat*diag(WK(:,:,k));
                    stimK=basisMat*(xK(:,k));

                    histEffect=exp(gamma*Hk')';
                    stimEffect=exp(stimK);
                    C = stimEffect.*histEffect;
                    M = 1./C;
                    lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                    ExpLambdaDelta = lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;
                    ExpLDSquaredTimesInvExp = (lambdaDelta).^2.*1./C;
                    ExpLDCubedTimesInvExpSquared = (lambdaDelta).^3.*M.^2 +Wk/2.*(3.*M.^4.*lambdaDelta.^3+12.*lambdaDelta.^3.*M.^3-12.*M.^4.*lambdaDelta.^4);

                    jacQ  = jacQ  - (Hk.*repmat(ExpLDSquaredTimesInvExp.*dN(k,:)',[1,size(Hk,2)]))'*Hk ...
                                  - (Hk.*repmat(ExpLDSquaredTimesInvExp,[1,size(Hk,2)]))'*Hk ...
                                  - (Hk.*repmat(2*ExpLDCubedTimesInvExpSquared,[1,size(Hk,2)]))'*Hk;

                 end


            end           

            Ic(1:R,1:R)=K*eye(size(Q))/(2*(Q)^2)+(eye(size(Q))/((Q)^3))*SumXkTerms;

            if(~isempty(windowTimes))
                Ic((R+1):(R+J),(R+1):(R+J)) = -jacQ;
            end
            xKDraw = zeros(numBasis,K,Mc);
            for r=1:numBasis  
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
            end



            Im=zeros(J+R,J+R);
            ImMC=zeros(J+R,J+R);

            for c=1:Mc

                gradQGammahat=zeros(size(gamma,2),1);
                gradQQhat=zeros(1,R);        
                if(strcmp(fitType,'poisson'))
                    for k=1:K
                        Hk=HkAll{k};
                        stimK=basisMat*(xKDraw(:,k,c));
                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK);
                        lambdaDelta = stimEffect.*histEffect;
                        gradQGammahat = gradQGammahat + Hk'*dN(k,:)' - Hk'*lambdaDelta;
                        if(k==1)
                            gradQQhat = ((xKDraw(:,k,c)-A*x0).*(xKDraw(:,k,c)-A*x0));
                        else
                            gradQQhat = gradQQhat+((xKDraw(:,k,c)-A*xKDraw(:,k-1,c)).*(xKDraw(:,k,c)-A*xKDraw(:,k-1,c)));
                        end

                    end
                elseif(strcmp(fitType,'binomial'))
                     for k=1:K
                        Hk=HkAll{k};
                        Wk = basisMat*diag(WK(:,:,k));
                        stimK=basisMat*(xKDraw(:,k,c));

                        histEffect=exp(gamma*Hk')';
                        stimEffect=exp(stimK);
    %                   
                        C = stimEffect.*histEffect;
                        M = 1./C;
                        lambdaDelta = exp(stimK+(gamma*Hk')')./(1+exp(stimK+(gamma*Hk')'));
                        ExpLambdaDelta = lambdaDelta+Wk.*(lambdaDelta.*(1-lambdaDelta).*(1-2*lambdaDelta))/2;
                        ExpLDSquaredTimesInvExp = (lambdaDelta).^2.*1./C;
                        ExpLDCubedTimesInvExpSquared = (lambdaDelta).^3.*M.^2 +Wk/2.*(3.*M.^4.*lambdaDelta.^3+12.*lambdaDelta.^3.*M.^3-12.*M.^4.*lambdaDelta.^4);


                        gradQGammahat = gradQGammahat + (Hk.*repmat(1-ExpLambdaDelta,[1,size(Hk,2)]))'*dN(k,:)' ...
                                          - (Hk.*repmat(ExpLDSquaredTimesInvExp./lambdaDelta,[1,size(Hk,2)]))'*lambdaDelta;
                        if(k==1)
                            gradQQhat = ((xKDraw(:,k,c)-A*x0).*(xKDraw(:,k,c)-A*x0));
                        else
                            gradQQhat = gradQQhat+((xKDraw(:,k,c)-A*xKDraw(:,k-1,c)).*(xKDraw(:,k,c)-A*xKDraw(:,k-1,c)));
                        end
                     end


                end

                gradQQhat = .5*eye(size(Q))/Q*gradQQhat - diag(K/2*eye(size(Q))/Q^2);
                ImMC(1:R,1:R)=ImMC(1:R,1:R)+gradQQhat*gradQQhat';
                if(~isempty(windowTimes))
                    ImMC((R+1):(R+J),(R+1):(R+J)) = ImMC((R+1):(R+J),(R+1):(R+J))+diag(diag(gradQGammahat*gradQGammahat'));
                end
            end
            Im=ImMC/Mc;

            InfoMatrix=Ic-Im; % Observed information matrix



        end
        function [spikeRateSig, ProbMat,sigMat]=computeSpikeRateCIs(xK,Wku,dN,t0,tf,fitType,delta,gamma,windowTimes,Mc,alphaVal)
             if(nargin<11 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<10 ||isempty(Mc))
                Mc=500;
            end

            [numBasis,K]=size(xK);

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end


    %         K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    Hk{k} = 0;
                end
                gamma=0;
            end

           for r=1:numBasis  
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
           end

           time=minTime:delta:maxTime;
           for c=1:Mc
               for k=1:K

                   if(strcmp(fitType,'poisson'))
                        stimK=basisMat*xKDraw(:,k,c);
                        histEffect=exp(gamma*Hk{k}')';
                        stimEffect=exp(stimK);
                        lambdaDelta(:,k,c) =stimEffect.*histEffect;
                   elseif(strcmp(fitType,'binomial'))
                        stimK=basisMat*xKDraw(:,k,c);
                        lambdaDelta(:,k,c)=exp(stimK+(gamma*Hk{k}')')./(1+exp(stimK+(gamma*Hk{k}')'));  
                   end  


               end
               lambdaC=Covariate(time,lambdaDelta(:,:,c)/delta,'\Lambda(t)');
               lambdaCInt= lambdaC.integral;
               spikeRate(c,:) = (1/(tf-t0))*(lambdaCInt.getValueAt(tf)-lambdaCInt.getValueAt(t0));

           end

           CIs = zeros(K,2);
           for k=1:K
               [f,x] = ecdf(spikeRate(:,k));
                CIs(k,1) = x(find(f<alphaVal,1,'last'));
                CIs(k,2) = x(find(f>(1-(alphaVal)),1,'first'));
           end
           spikeRateSig = Covariate(1:K, mean(spikeRate),['(' num2str(tf) '-' num2str(t0) ')^-1 * \Lambda(' num2str(tf) '-' num2str(t0) ')'],'Trial','k','Hz');
           ciSpikeRate = ConfidenceInterval(1:K,CIs,'CI_{spikeRate}','Trial','k','Hz');
           spikeRateSig.setConfInterval(ciSpikeRate);


           if(nargout>1)
               ProbMat = zeros(K,K);
               for k=1:K
                   for m=(k+1):K

                       ProbMat(k,m)=sum(spikeRate(:,m)>spikeRate(:,k))./Mc;
                   end
               end
           end


           if(nargout>2)
                sigMat= double(ProbMat>(1-alphaVal));
           end


        end
        function [spikeRateSig, ProbMat,sigMat]=computeSpikeRateDiffCIs(xK,Wku,dN,time1,time2,fitType,delta,gamma,windowTimes,Mc,alphaVal)
             if(nargin<11 ||isempty(alphaVal))
                alphaVal =.05;
            end
            if(nargin<10 ||isempty(Mc))
                Mc=500;
            end

            [numBasis,K]=size(xK);

            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            if(~isempty(numBasis))
                basisWidth = (maxTime-minTime)/numBasis;
                sampleRate=1/delta;
                unitPulseBasis=nstColl.generateUnitImpulseBasis(basisWidth,minTime,maxTime,sampleRate);
                basisMat = unitPulseBasis.data;
            end


    %         K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
                    Hk{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
                    Hk{k} = 0;
                end
                gamma=0;
            end

           for r=1:numBasis  
                WkuTemp=squeeze(Wku(r,r,:,:));
    %             [vec,val]=eig(Wku ); val(val<=0)=eps;
    %             Wku =vec*val*vec';
                [chol_m,p]=chol(WkuTemp);
                if(numel(chol_m)==1)
                    chol_m = diag(repmat(chol_m,[K 1]));
                end
                for c=1:Mc % for r-th step function simulate the path of size K
                    z=zeros(K,1);
                    z=normrnd(0,1,K,1);
                    xKDraw(r,:,c)=xK(r,:)+(chol_m'*z)';
                end
           end

           timeWindow=minTime:delta:maxTime;
           for c=1:Mc
               for k=1:K

                   if(strcmp(fitType,'poisson'))
                        stimK=basisMat*xKDraw(:,k,c);
                        histEffect=exp(gamma*Hk{k}')';
                        stimEffect=exp(stimK);
                        lambdaDelta(:,k,c) =stimEffect.*histEffect;
                   elseif(strcmp(fitType,'binomial'))
                        stimK=basisMat*xKDraw(:,k,c);
                        lambdaDelta(:,k,c)=exp(stimK+(gamma*Hk{k}')')./(1+exp(stimK+(gamma*Hk{k}')'));  
                   end  


               end
               lambdaC=Covariate(timeWindow,lambdaDelta(:,:,c)/delta,'\Lambda(t)');
               lambdaCInt= lambdaC.integral;
               spikeRate(c,:) = (1/(max(time1)-min(time1)))*(lambdaCInt.getValueAt(max(time1))-lambdaCInt.getValueAt(min(time1))) ...
                                - (1/(max(time2)-min(time2)))*(lambdaCInt.getValueAt(max(time2))-lambdaCInt.getValueAt(min(time2)));


           end

           CIs = zeros(K,2);
           for k=1:K
               [f,x] = ecdf(spikeRate(:,k));
                CIs(k,1) = x(find(f<alphaVal,1,'last')); %not alpha/2 since this is a once sided comparison
                CIs(k,2) = x(find(f>(1-(alphaVal)),1,'first'));
           end
           spikeRateSig = Covariate(1:K, mean(spikeRate),['(t_{1f}-t_{1o})^-1 * \Lambda(t_{1f}-t_{1o}) - (t_{2f}-t_{2o})^-1 * \Lambda(t_{2f}-t_{2o}) '],'Trial','k','Hz');
           ciSpikeRate = ConfidenceInterval(1:K,CIs,'CI_{spikeRate}','Trial','k','Hz');
           spikeRateSig.setConfInterval(ciSpikeRate);


           if(nargout>1)
               ProbMat = zeros(K,K);
               for k=1:K
                   for m=(k+1):K

                       ProbMat(k,m)=sum(spikeRate(:,m)>spikeRate(:,k))./Mc;
                   end
               end
           end


           if(nargout>2)
                sigMat= double(ProbMat>(1-alphaVal));
           end


        end
        
        %% Kalman Filter EM
        function C = KF_EMCreateConstraints(EstimateA, AhatDiag,QhatDiag,QhatIsotropic,RhatDiag,RhatIsotropic,Estimatex0,EstimatePx0, Px0Isotropic,mcIter,EnableIkeda)
            %By default, all parameters are estimated. To empose diagonal
            %structure on the EM parameter results must pass in the
            %constraints element
            if(nargin<11 || isempty(EnableIkeda))
                EnableIkeda=0;
            end
            if(nargin<10 || isempty(mcIter))
                mcIter=1000;
            end
            if(nargin<9 || isempty(Px0Isotropic))
                Px0Isotropic=0;
            end
            if(nargin<8 || isempty(EstimatePx0))
                EstimatePx0=1;
            end
            if(nargin<7 || isempty(Estimatex0))
                Estimatex0=1;
            end
            if(nargin<6 || isempty(RhatIsotropic))
                RhatIsotropic=0;
            end
            if(nargin<5 || isempty(RhatDiag))
                RhatDiag=1;
            end
            if(nargin<4 || isempty(QhatIsotropic))
                QhatIsotropic=0;
            end
            if(nargin<3 || isempty(QhatDiag))
                QhatDiag=1;
            end
            if(nargin<2)
                AhatDiag=0;
            end
            if(nargin<1)
                EstimateA=1;
            end
            C.EstimateA = EstimateA;
            C.AhatDiag = AhatDiag;
            C.QhatDiag = QhatDiag;
            if(QhatDiag && QhatIsotropic)
                C.QhatIsotropic=1;
            else
                C.QhatIsotropic=0;
            end
            C.RhatDiag = RhatDiag;
            if(RhatDiag && RhatIsotropic)
                C.RhatIsotropic=1;
            else
                C.RhatIsotropic=0;
            end
            C.Estimatex0 = Estimatex0;
            C.EstimatePx0 = EstimatePx0;
            if(EstimatePx0 && Px0Isotropic)
                C.Px0Isotropic=1;
            else
                C.Px0Isotropic=0; 
            end
            C.mcIter = mcIter;
            C.EnableIkeda = EnableIkeda;
        end
        function [xKFinal,WKFinal,Ahat, Qhat, Chat, Rhat,alphahat, x0hat, Px0hat, IC, SE, Pvals, nIter]=KF_EM(y, Ahat0, Qhat0, Chat0, Rhat0, alphahat0, x0, Px0,KFEM_Constraints)
            numStates = size(Ahat0,1);
            
            if(nargin<9 || isempty(KFEM_Constraints))
                KFEM_Constraints=DecodingAlgorithms.KF_EMCreateConstraints;
            end
            if(nargin<8 || isempty(Px0))
                Px0=10e-10*eye(numStates,numStates);
            end
            if(nargin<7 || isempty(x0))
                x0=zeros(numStates,1);
            end
            
    %         tol = 1e-3; %absolute change;
            tolAbs = 1e-3;
            tolRel = 1e-3;
            llTol  = 1e-3;
            cnt=1;

            maxIter = 100;

            
            A0 = Ahat0;
            Q0 = Qhat0;
            C0 = Chat0;
            R0 = Rhat0;
            alpha0 = alphahat0;
           
           
            Ahat{1} = A0;
            Qhat{1} = Q0;
            Chat{1} = C0;
            Rhat{1} = R0;
            x0hat{1} = x0;
            Px0hat{1} = Px0;
            alphahat{1} = alpha0;
            yOrig=y;
            numToKeep=10;
            scaledSystem=1;
            
            if(scaledSystem==1)
                Tq = eye(size(Qhat{1}))/(chol(Qhat{1}));
                Tr = eye(size(Rhat{1}))/(chol(Rhat{1}));
                Ahat{1}= Tq*Ahat{1}/Tq;
                Chat{1}= Tr*Chat{1}/Tq;
                Qhat{1}= Tq*Qhat{1}*Tq';
                Rhat{1}= Tr*Rhat{1}*Tr';
                y= Tr*y;
                x0hat{1} = Tq*x0;
                Px0hat{1} = Tq*Px0*Tq';
                alphahat{1}= Tr*alphahat{1};  
            end

            cnt=1;
            dLikelihood(1)=inf;
            negLL=0;
            IkedaAcc=KFEM_Constraints.EnableIkeda;
            %Forward EM
            stoppingCriteria =0;

            disp('                       Kalman Filter/Gaussian Observation EM Algorithm                        ');     
            while(stoppingCriteria~=1 && cnt<=maxIter)
                 storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
                 storeIndP1= mod(cnt,numToKeep)+1;
                 storeIndM1= mod(cnt-2,numToKeep)+1;
                disp('--------------------------------------------------------------------------------------------------------');
                disp(['Iteration #' num2str(cnt)]);
                disp('--------------------------------------------------------------------------------------------------------');
                
                
                [x_K{storeInd},W_K{storeInd},ll(cnt),ExpectationSums{storeInd}]=...
                    DecodingAlgorithms.KF_EStep(Ahat{storeInd},Qhat{storeInd},Chat{storeInd},Rhat{storeInd}, y, alphahat{storeInd}, x0hat{storeInd}, Px0hat{storeInd});
                
                [Ahat{storeIndP1}, Qhat{storeIndP1}, Chat{storeIndP1}, Rhat{storeIndP1}, alphahat{storeIndP1},x0hat{storeIndP1},Px0hat{storeIndP1}] ...
                    = DecodingAlgorithms.KF_MStep(y,x_K{storeInd},x0hat{storeInd}, Px0hat{storeInd},ExpectationSums{storeInd},KFEM_Constraints);
              
                if(IkedaAcc==1)
                    disp(['****Ikeda Acceleration Step****']);
                    %y=Cx+alpha+wk wk~Normal with covariance Rk
                     ykNew = mvnrnd((Chat{storeIndP1}*x_K{storeInd}+alphahat{storeIndP1}*ones(1,size(x_K{storeInd},2)))',Rhat{storeIndP1})';
                     
                                    
                     [x_KNew,W_KNew,llNew,ExpectationSumsNew]=...
                        DecodingAlgorithms.KF_EStep(Ahat{storeInd},Qhat{storeInd},Chat{storeInd},Rhat{storeInd}, ykNew, alphahat{storeInd},x0, Px0);
         
                     [AhatNew, QhatNew, ChatNew, RhatNew, alphahatNew,x0new,Px0new] ...
                        = DecodingAlgorithms.KF_MStep(ykNew,x_KNew, x0hat{storeInd}, Px0hat{storeInd}, ExpectationSumsNew,KFEM_Constraints);
               
                    Ahat{storeIndP1} = 2*Ahat{storeIndP1}-AhatNew;
                    Qhat{storeIndP1} = 2*Qhat{storeIndP1}-QhatNew;
                    Qhat{storeIndP1} = (Qhat{storeIndP1}+Qhat{storeIndP1}')/2;
                    Chat{storeIndP1} = 2*Chat{storeIndP1}-ChatNew;
                    Rhat{storeIndP1} = 2*Rhat{storeIndP1}-RhatNew;
                    Rhat{storeIndP1} = (Rhat{storeIndP1}+Rhat{storeIndP1}')/2;
                    alphahat{storeIndP1}=2*alphahat{storeIndP1}-alphahatNew;
                    
%                     x0hat{storeIndP1}   = 2*x0hat{storeIndP1} - x0new;
%                     Px0hat{storeIndP1}  = 2*Px0hat{storeIndP1}- Px0new;
%                     [V,D] = eig(Px0hat{storeIndP1});
%                     D(D<0)=1e-9;
%                     Px0hat{storeIndP1} = V*D*V';
%                     Px0hat{storeIndP1}  = (Px0hat{storeIndP1}+Px0hat{storeIndP1}')/2;
                    
               
                end
                if(KFEM_Constraints.EstimateA==0)
                    Ahat{storeIndP1}=Ahat{storeInd};
                end
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(ll(cnt)-ll(cnt-1));%./abs(logll(cnt-1));
                end
                if(cnt==1)
                    QhatInit = Qhat{1};
                    RhatInit = Rhat{1};
                    xKInit = x_K{1};
                end
                %Plot the progress
%                 if(mod(cnt,2)==0)
                if(cnt==1)
                    scrsz = get(0,'ScreenSize');
                    h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
                end
                    figure(h);
                    time = 0:(size(y,2)-1);
                    
                    subplot(2,5,[1 2 6 7]); plot(1:cnt,ll,'k','Linewidth', 2); hy=ylabel('Log Likelihood'); hx=xlabel('Iteration'); axis auto;
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    subplot(2,5,3:5); hNew=plot(time, x_K{storeInd}','Linewidth', 2); hy=ylabel('States'); hx=xlabel('time [s]');
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(time, xKInit','--','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    subplot(2,5,8); hNew=plot(diag(Qhat{storeInd}),'o','Linewidth', 2); hy=ylabel('Q'); hx=xlabel('Diagonal Entry');
                    set(gca, 'XTick'       , 1:1:length(diag(Qhat{storeInd})));
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(diag(QhatInit),'r.','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    subplot(2,5,9); hNew=plot(diag(Rhat{storeInd}),'o','Linewidth', 2); hy=ylabel('R'); hx=xlabel('Diagonal Entry');
                    set(gca, 'XTick'       , 1:1:length(diag(Rhat{storeInd})));
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(diag(RhatInit),'r.','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    
                    subplot(2,5,10); imagesc(Rhat{storeInd}); ht=title('R Matrix Image'); 
                    set(gca, 'XTick'       , 1:1:length(diag(Rhat{storeInd})), 'YTick', 1:1:length(diag(Rhat{storeInd})));
                    set(ht,'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    drawnow;
                    hold off;
%                 end
                
                if(cnt==1)
                    dMax=inf;
                else
                 dQvals = max(max(abs(sqrt(Qhat{storeInd})-sqrt(Qhat{storeIndM1}))));
                 dRvals = max(max(abs(sqrt(Rhat{storeInd})-sqrt(Rhat{storeIndM1}))));
                 dAvals = max(max(abs((Ahat{storeInd})-(Ahat{storeIndM1}))));
                 dCvals = max(max(abs((Chat{storeInd})-(Chat{storeIndM1}))));
                 dAlphavals = max(abs((alphahat{storeInd})-(alphahat{storeIndM1})));
                 dMax = max([dQvals,dRvals,dAvals,dCvals,dAlphavals]);
                end

% 
%                 dQRel = max(abs(dQvals./sqrt(Qhat(:,storeIndM1))));
%                 dGammaRel = max(abs(dGamma./gammahat(storeIndM1,:)));
%                 dMaxRel = max([dQRel,dGammaRel]);

                if(cnt==1)
                    disp(['Max Parameter Change: N/A']);
                else
                    disp(['Max Parameter Change: ' num2str(dMax)]);
                end 
                cnt=(cnt+1);
                
                if(dMax<tolAbs)
                    stoppingCriteria=1;
                    display(['         EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
                    negLL=0;
                end
            
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
                    stoppingCriteria=1;
                    display(['         EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
                    negLL=1;
                end
            
                
            end
            disp('--------------------------------------------------------------------------------------------------------');
            
            maxLLIndex  = find(ll == max(ll),1,'first');
            maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
            if(maxLLIndex==1)
%                 maxLLIndex=cnt-1;
                maxLLIndex =1;
                maxLLIndMod = 1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
               maxLLIndMod = 1;
%             else
%                maxLLIndMod = mod(maxLLIndex,numToKeep); 
               
            end
            nIter   = cnt-1;  
%             maxLLIndMod
           
            xKFinal = x_K{maxLLIndMod};
            WKFinal = W_K{maxLLIndMod};
            Ahat = Ahat{maxLLIndMod};
            Qhat = Qhat{maxLLIndMod};
            Chat = Chat{maxLLIndMod};
            Rhat = Rhat{maxLLIndMod};
            alphahat = alphahat{maxLLIndMod};
            x0hat =x0hat{maxLLIndMod};
            Px0hat=Px0hat{maxLLIndMod};
            
             if(scaledSystem==1)
               Tq = eye(size(Qhat))/(chol(Q0));
               Tr = eye(size(Rhat))/(chol(R0));
               Ahat=Tq\Ahat*Tq;
               Qhat=(Tq\Qhat)/Tq';
               Chat=Tr\Chat*Tq;
               Rhat=(Tr\Rhat)/Tr';
               alphahat=Tr\alphahat;
               xKFinal = Tq\xKFinal;
               x0hat = Tq\x0hat;
               Px0hat= (Tq\Px0hat)/(Tq');
               tempWK =zeros(size(WKFinal));
               for kk=1:size(WKFinal,3)
                tempWK(:,:,kk)=(Tq\WKFinal(:,:,kk))/Tq';
               end
               WKFinal = tempWK;
             end
            
            ll = ll(maxLLIndex);
            ExpectationSumsFinal = ExpectationSums{maxLLIndMod};

            if(nargout>10)
                [SE, Pvals]=DecodingAlgorithms.KF_ComputeParamStandardErrors(y, xKFinal, WKFinal, Ahat, Qhat, Chat, Rhat, alphahat, x0hat, Px0hat, ExpectationSumsFinal, KFEM_Constraints);
            end
             %Compute number of parameters
            if(KFEM_Constraints.EstimateA==1 && KFEM_Constraints.AhatDiag==1)
                n1=size(Ahat,1); 
            elseif(KFEM_Constraints.EstimateA==1 && KFEM_Constraints.AhatDiag==0)
                n1=numel(Ahat);
            else 
                n1=0;
            end
            if(KFEM_Constraints.QhatDiag==1 && KFEM_Constraints.QhatIsotropic==1)
                n2=1;
            elseif(KFEM_Constraints.QhatDiag==1 && KFEM_Constraints.QhatIsotropic==0)
                n2=size(Qhat,1);
            else
                n2=numel(Qhat);
            end

            n3=numel(Chat); 
            if(KFEM_Constraints.RhatDiag==1 && KFEM_Constraints.RhatIsotropic==1)
                n4=1;
            elseif(KFEM_Constraints.QhatDiag==1 && KFEM_Constraints.QhatIsotropic==0)
                n4=size(Rhat,1);
            else
                n4=numel(Rhat);
            end

            if(KFEM_Constraints.EstimatePx0==1 && KFEM_Constraints.Px0Isotropic==1)
                n5=1;
            elseif(KFEM_Constraints.EstimatePx0==1 && KFEM_Constraints.Px0Isotropic==0)
                n5=size(Px0hat,1);
            else
                n5=0;
            end

            if(KFEM_Constraints.Estimatex0==1)   
                n6=size(x0hat,1);
            else
                n6=0;
            end

            n7=size(alphahat,1);
           
            nTerms=n1+n2+n3+n4+n5+n6+n7;
            K  = size(y,2); 
            Dx = size(Ahat,2);
            sumXkTerms = ExpectationSums{maxLLIndMod}.sumXkTerms;
            llobs = ll + Dx*K/2*log(2*pi)+K/2*log(det(Qhat))...
                + 1/2*trace(Qhat\sumXkTerms)...
                + Dx/2*log(2*pi)+1/2*log(det(Px0hat)) ...
                + 1/2*Dx;
            AIC = 2*nTerms - 2*llobs;
            AICc= AIC+ 2*nTerms*(nTerms+1)/(K-nTerms-1);
            BIC = -2*llobs+nTerms*log(K);
            IC.AIC = AIC;
            IC.AICc= AICc;
            IC.BIC = BIC;
            IC.llobs = llobs;
            IC.llcomp=ll;
            
            
            
        end
        function [SE,Pvals,nTerms] = KF_ComputeParamStandardErrors(y, xKFinal, WKFinal, Ahat, Qhat, Chat, Rhat, alphahat, x0hat, Px0hat, ExpectationSumsFinal, KFEM_Constraints)
         % Use inverse observed information matrix to estimate the standard errors of the estimated model parameters
         % Requires computation of the complete information matrix and an estimate of the missing information matrix

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Complete Information Matrices     
        % Recall from McLachlan and Krishnan Eq. 4.7
        %    Io(theta;y) = Ic(theta;y) - Im(theta;y)
        %    Io(theta;y) = Ic(theta;y) - cov(Sc(X;theta)Sc(X;theta)')
        % where Sc(X;theta) is the score vector of the complete log likelihood
        % function evaluated at theta. We first compute Ic term by term and then
        % approximate the covariance term using Monte Carlo approximation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(nargin<12 || isempty(KFEM_Constraints))
                KFEM_Constraints=DecodingAlgorithms.KF_EMCreateConstraints;
            end

            if(KFEM_Constraints.AhatDiag==1)
                IAComp=zeros(numel(diag(Ahat)),numel(diag(Ahat)));
            else
                IAComp=zeros(numel(Ahat),numel(Ahat));
            end
            [n1,n2] =size(Ahat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            N=size(y,2);

            if(KFEM_Constraints.AhatDiag==1)
                for l=1:n1
                    for m=l
                        termMat=Qhat\el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1.*eye(n1,n2);
                        termvec = diag(termMat);
                        IAComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            else
                for l=1:n1
                    for m=1:n2
                        termMat=(inv(Qhat))*el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1;
                        termvec=reshape(termMat',1,numel(Ahat));
                        IAComp(:,cnt)=termvec';
                        cnt=cnt+1;
                    end
                end
            end


            ICComp=zeros(numel(Chat),numel(Chat));
            [n1,n2] =size(Chat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            for l=1:n1
                for m=1:n2
                    termMat=Rhat\el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkxk;
                    termvec=reshape(termMat',1,numel(Chat));
                    ICComp(:,cnt)=termvec';
                    cnt=cnt+1;
                end
            end

            [n1,n2] =size(Rhat);
            ei=(eye(n1,n1));
            ej=(eye(n2,n2));
            cnt=1;
            [dy,N]=size(y);
            dx=size(xKFinal,1);

        %     if(KFEM_Constraints.RhatDiag==1)
        %         if(KFEM_Constraints.RhatIsotropic)
        %             IRinvComp=zeros(1,1);
        %             IRinvComp =  0.5*N*dy*Rhat(1,1)^2; 
        %         else
        %             IRinvComp=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
        %             for i=1:n1
        %                 for j=i
        %                     termMat= N/2*(Rhat)\ei(:,i)*ej(:,j)'/(Rhat);
        %                     termvec=diag(termMat);
        %         %             termvec=reshape(termMat',1,numel(Rhat));
        %                     IRinvComp(cnt,:)=termvec';
        %                     cnt=cnt+1;
        %                 end
        %             end
        %         end
        %     else
        %         IRinvComp=zeros(numel(Rhat),numel(Rhat));
        %         for i=1:n1
        %             for j=1:n2
        %                 termMat= N/2*(Rhat)\ei(:,i)*ej(:,j)'/(Rhat);
        %                 termvec=reshape(termMat',1,numel(Rhat));
        %                 IRinvComp(cnt,:)=termvec;
        %                 cnt=cnt+1;
        %             end
        %         end
        %     end

            [n1,n2] =size(Rhat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            N=size(y,2);
            if(KFEM_Constraints.RhatDiag==1)
                if(KFEM_Constraints.RhatIsotropic==1)
                    IRComp = 0.5*N*dy*Rhat(1,1)^(-2);
                else
                    IRComp=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
                    for l=1:n1
                        for m=l
                            termMat= N/2*(Rhat)\em(:,m)*el(:,l)'/(Rhat);
                            termvec=diag(termMat);
                            IRComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            else
                IRComp=zeros(numel((Rhat)),numel((Rhat)));
                for l=1:n1
                    for m=1:n2
                        termMat= N/2*(Rhat)\em(:,m)*el(:,l)'/(Rhat);
                        termvec=reshape(termMat',1,numel(Rhat));
                        IRComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            end

        %     [n1,n2] =size(Qhat);
        %     ei=(eye(n1,n1));
        %     ej=(eye(n2,n2));
        %     cnt=1;
        %     N=size(y,2);
        %     dx=size(xKFinal,1);
        %     if(KFEM_Constraints.QhatDiag==1)
        %         if(KFEM_Constraints.QhatIsotropic==1)
        %             IQinvComp=zeros(1,1);
        %             IQinvComp =  0.5*N*dx*Qhat(1,1)^2; 
        %             
        %         else
        %             IQinvComp=zeros(numel(diag(Qhat)),numel(diag(Qhat)));
        %             for i=1:n1
        %                 for j=i
        %                     termMat= N/2*(Qhat)\ei(:,i)*ej(:,j)'/(Qhat);
        %                     termvec=diag(termMat);
        %         %             termvec=reshape(termMat',1,numel(Qhat));
        %                     IQinvComp(cnt,:)=termvec';
        %                     cnt=cnt+1;
        %                 end
        %             end
        %         end
        %     else
        %         IQinvComp=zeros(numel(Qhat),numel(Qhat));
        %         for i=1:n1
        %             for j=1:n2
        %                 termMat= N/2*(Qhat)\ei(:,i)*ej(:,j)'/(Qhat);
        %                 termvec=reshape(termMat',1,numel(Qhat));
        %                 IQinvComp(cnt,:)=termvec';
        %                 cnt=cnt+1;
        %             end
        %         end
        %     end

            [n1,n2] =size(Qhat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            N=size(y,2);
            if(KFEM_Constraints.QhatDiag==1)
                if(KFEM_Constraints.QhatIsotropic==1)
                    IQComp=zeros(1,1);
                    IQComp =  0.5*N*dx*Qhat(1,1)^(-2); 
                else
                    IQComp=zeros(numel(diag(Qhat)),numel(diag(Qhat)));
                    for l=1:n1
                        for m=l
                            termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                            termvec=diag(termMat);
                            IQComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            else
                IQComp=zeros(numel(Qhat),numel(Qhat));
                for l=1:n1
                    for m=1:n2
                        termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                        termvec=reshape(termMat',1,numel(Qhat));
                        IQComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            end

            if(KFEM_Constraints.EstimatePx0==1)
                if(KFEM_Constraints.Px0Isotropic==1)
        %             ISinvComp =  0.5*dx*Px0hat(1,1)^2;
                    ISComp =  0.5*dx*Px0hat(1,1)^(-2);
                else
        %             ISinvComp=zeros(numel(diag(Px0hat)),numel(diag(Px0hat)));
        %             [n1,n2] =size(Px0hat);
        %             ei=(eye(n1,n1));
        %             ej=(eye(n2,n2));
        %             cnt=1;
        %             for i=1:n1
        %                 for j=i
        %                     termMat= 1/2*(Px0hat)\ei(:,i)*ej(:,j)'/(Px0hat);
        %         %             termvec=reshape(termMat',1,numel(Px0hat));
        %                     termvec=diag(termMat);
        %                     ISinvComp(cnt,:)=termvec';
        %                     cnt=cnt+1;
        %                 end
        %             end

                    ISComp=zeros(numel(diag(Px0hat)),numel(diag(Px0hat)));
                    [n1,n2] =size(Px0hat);
                    el=(eye(n1,n1));
                    em=(eye(n2,n2));
                    cnt=1;
                    for l=1:n1
                        for m=l
                            termMat= 1/2*(Px0hat)\em(:,m)*el(:,l)'/(Px0hat);
                            termvec=diag(termMat);
                %             termvec=reshape(termMat',1,numel(Rhat));
                            ISComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            end

            if(KFEM_Constraints.Estimatex0==1)
                Ix0Comp=eye(size(Px0hat))/Px0hat+(Ahat'/Qhat)*Ahat;
            end

            IAlphaComp = N*eye(size(Rhat))/Rhat;
            if(KFEM_Constraints.EstimateA==1)
                n1=size(IAComp,1); 
            else
                n1=0;
            end
            n2=size(IQComp,1); n3=size(ICComp,1);
            n4=size(IRComp,1); 
            if(KFEM_Constraints.EstimatePx0==1)
                n5=size(ISComp,1); 
            else
                n5=0;
            end
            if(KFEM_Constraints.Estimatex0==1)   
                n6=size(Ix0Comp,1);
            else
                n6=0;
            end
            n7=size(IAlphaComp,1);
            nTerms=n1+n2+n3+n4+n5+n6+n7;
            IComp = zeros(nTerms,nTerms);
            if(KFEM_Constraints.EstimateA==1)
                IComp(1:n1,1:n1)=IAComp;
            end
            offset=n1+1;
            IComp(offset:(n1+n2),offset:(n1+n2))=IQComp;
            offset=n1+n2+1;
            IComp(offset:(n1+n2+n3),offset:(n1+n2+n3))=ICComp;
            offset=n1+n2+n3+1;
            IComp(offset:(n1+n2+n3+n4),offset:(n1+n2+n3+n4))=IRComp;
            offset=n1+n2+n3+n4+1;
            if(KFEM_Constraints.EstimatePx0==1);
                IComp(offset:(n1+n2+n3+n4+n5),offset:(n1+n2+n3+n4+n5))=ISComp;
            end
            offset=n1+n2+n3+n4+n5+1;
            if(KFEM_Constraints.Estimatex0==1)
                IComp(offset:(n1+n2+n3+n4+n5+n6),offset:(n1+n2+n3+n4+n5+n6))=Ix0Comp;
            end
            offset=n1+n2+n3+n4+n5+n6+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7),offset:(n1+n2+n3+n4+n5+n6+n7))=IAlphaComp;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Missing Information Matrix
            %Approximate cov(Sc(X;theta)Sc(X;theta)')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K=size(y,2);
            Mc=KFEM_Constraints.mcIter;
            xKDraw = zeros(size(xKFinal,1),N,Mc);
        %     IMissing = zeros(nTerms,nTerms);
        % 
        %     if(KFEM_Constraints.AhatDiag==1)
        %         ScoreAMc = zeros(numel(diag(Ahat)),Mc);
        %         IAMissing=zeros(numel(diag(Ahat)),numel(diag(Ahat)));
        %     else
        %         ScoreAMc = zeros(numel(Ahat),Mc);
        %         IAMissing=zeros(numel(Ahat),numel(Ahat));
        %     end
        % 
        %     ScoreCMc = zeros(numel(Chat),Mc);
        %     if(KFEM_Constraints.RhatDiag==1)
        %         if(KFEM_Constraints.RhatIsotropic==1)
        %             ScoreRinvMc = zeros(1,Mc);  
        %             ScoreRMc = zeros(1,Mc);
        %             IRinvMissing=zeros(1,1);
        %             IRMissing=zeros(1,1);
        %         else
        %             ScoreRinvMc = zeros(numel(diag(Rhat)),Mc);  
        %             ScoreRMc = zeros(numel(diag(Rhat)),Mc);
        %             IRinvMissing=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
        %             IRMissing=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
        %         end
        %     else
        %         ScoreRMc = zeros(numel(Rhat),Mc);
        %         ScoreRinvMc = zeros(numel(Rhat),Mc);
        %         IRMissing=zeros(numel(Rhat),numel(Rhat));
        %         IRinvMissing=zeros(numel(Rhat),numel(Rhat));
        %     end
        % 
        %     if(KFEM_Constraints.QhatDiag==1)
        %         if(KFEM_Constraints.QhatIsotropic==1)
        %             ScoreQinvMc = zeros(1,Mc);
        %             ScoreQMc    = zeros(1,Mc); 
        %         else
        %             ScoreQinvMc = zeros(numel(diag(Qhat)),Mc);
        %             ScoreQMc    = zeros(numel(diag(Qhat)),Mc);
        %         end
        %     else
        %         ScoreQMc = zeros(numel(Qhat),Mc);
        %     end
        %     ScoreAlphaMc = zeros(numel(alphahat),Mc);

            % Generate the Monte Carlo samples for the unobserved data
            for n=1:N
                WuTemp=(WKFinal(:,:,n));
                [chol_m,p]=chol(WuTemp);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                xKDraw(:,n,:)=repmat(xKFinal(:,n),[1 Mc])+(chol_m*z);
            end


            if(KFEM_Constraints.EstimatePx0|| KFEM_Constraints.Estimatex0)
                [chol_m,p]=chol(Px0hat);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                x0Draw=repmat(x0hat,[1 Mc])+(chol_m*z); 
            else
               x0Draw=repmat(x0hat, [1 Mc]);

            end

            IMc = zeros(nTerms,nTerms,Mc);
            % Emperically estimate the covariance of the score
        %     if matlabpool('size') == 0 % checking to see if my pool is already open
        %         matlabpool;
        %     end
            pools = matlabpool('size'); %number of parallel workers 
            if(pools==0) % parallel toolbox is not enabled;
                for c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Dy=size(y,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);
                    Sykyk = zeros(Dy,Dy);
                    Sxkyk = zeros(Dx,Dy);        

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                        Sykyk = Sykyk+(y(:,k)-alphahat)*(y(:,k)-alphahat)';
                        Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alphahat)';

                    end
                    Sx0x0 = x_0*x_0';
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    Sykyk = 0.5*(Sykyk+Sykyk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    sumYkTerms = Sykyk - Chat*Sxkyk - Sxkyk'*Chat' + Chat*Sxkxk*Chat';      
                    Sxkxkm1 = Sxkm1xk';
                    Sykxk = Sxkyk';


                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    sumYkTerms=0.5*(sumYkTerms+sumYkTerms');
                    if(KFEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(KFEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    end

                    ScorC=Rhat\(Sykxk-Chat*Sxkxk);
                    ScoreCMc=reshape(ScorC',numel(ScorC),1);

                    if(KFEM_Constraints.QhatDiag)
                        if(KFEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end


                    ScoreAlphaMc = sum(Rhat\(y-Chat*x_K-alphahat*ones(1,N)),2);
                    if(KFEM_Constraints.RhatDiag)
                        if(KFEM_Constraints.RhatIsotropic)
                            ScoreR  =-.5*(K*Dy*Rhat(1,1)^(-1) - Rhat(1,1)^(-2)*trace(sumYkTerms));
                        else
                            ScoreR  =(-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat)));
                        end
                        ScoreRMc = diag(ScoreR);
                    else
                        ScoreR   =-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat));
                        ScoreRMc =reshape(ScoreR',numel(ScoreR),1);
                    end


                    if(KFEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);

                    if(KFEM_Constraints.EstimateA==1)
                        ScoreVec = ScoreAMc;
                    else
                        ScoreVec = [];
                    end
                    ScoreVec = [ScoreVec; ScoreQMc; ScoreCMc; ScoreRMc];
                    if(KFEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(KFEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreAlphaMc];
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end
            else %Use the parallel toolbox
                parfor c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Dy=size(y,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);
                    Sykyk = zeros(Dy,Dy);
                    Sxkyk = zeros(Dx,Dy);        

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                        Sykyk = Sykyk+(y(:,k)-alphahat)*(y(:,k)-alphahat)';
                        Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alphahat)';

                    end
                    Sx0x0 = x_0*x_0';
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    Sykyk = 0.5*(Sykyk+Sykyk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    sumYkTerms = Sykyk - Chat*Sxkyk - Sxkyk'*Chat' + Chat*Sxkxk*Chat';      
                    Sxkxkm1 = Sxkm1xk';
                    Sykxk = Sxkyk';


                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    sumYkTerms=0.5*(sumYkTerms+sumYkTerms');
                    if(KFEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(KFEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    else 
                        ScoreAMc=[];
                    end

                    ScorC=Rhat\(Sykxk-Chat*Sxkxk);
                    ScoreCMc=reshape(ScorC',numel(ScorC),1);

                    if(KFEM_Constraints.QhatDiag)
                        if(KFEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end


                    ScoreAlphaMc = sum(Rhat\(y-Chat*x_K-alphahat*ones(1,N)),2);
                    if(KFEM_Constraints.RhatDiag)
                        if(KFEM_Constraints.RhatIsotropic)
                            ScoreR  =-.5*(K*Dy*Rhat(1,1)^(-1) - Rhat(1,1)^(-2)*trace(sumYkTerms));
                        else
                            ScoreR  =(-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat)));
                        end
                        ScoreRMc = diag(ScoreR);
                    else
                        ScoreR   =-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat));
                        ScoreRMc =reshape(ScoreR',numel(ScoreR),1);
                    end


                    if(KFEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);
                    ScoreVec = [ScoreAMc; ScoreQMc; ScoreCMc; ScoreRMc];
                    if(KFEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(KFEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreAlphaMc];
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end

            end
            IMissing = 1/Mc*sum(IMc,3);
            IObs  = IComp-IMissing;  
            invIObs = eye(size(IObs))/IObs;
        %     figure(1); subplot(1,2,1); imagesc(invIObs); subplot(1,2,2); imagesc(nearestSPD(invIObs));
            invIObs = nearestSPD(invIObs); % Find the nearest positive semidefinite approximation for the variance matrix
            VarVec = (diag(invIObs));
            SEVec = sqrt(VarVec);
            SEAterms = SEVec(1:n1);
            SEQterms = SEVec(n1+1:(n1+n2));
            SECterms = SEVec(n1+n2+1:(n1+n2+n3));
            SERterms = SEVec(n1+n2+n3+1:(n1+n2+n3+n4));
            SEPx0terms=SEVec(n1+n2+n3+n4+1:(n1+n2+n3+n4+n5));
            SEx0terms=SEVec(n1+n2+n3+n4+n5+1:(n1+n2+n3+n4+n5+n6));
            SEAlphaterms=SEVec(n1+n2+n3+n4+n5+n6+1:(n1+n2+n3+n4+n5+n6+n7));

        %     matlabpool close;

        %     figure(1); 
        %     subplot(1,3,1); image(IObs); subplot(1,3,2); image(IComp); subplot(1,3,3); image(IMissing);
            if(KFEM_Constraints.EstimatePx0==1)
                SES = diag(SEPx0terms);
            end
            if(KFEM_Constraints.Estimatex0==1)
                SEx0=SEx0terms;
            end
            if(KFEM_Constraints.EstimateA==1)
                if(KFEM_Constraints.AhatDiag==1)
                    SEA=diag(SEAterms);
                else
                    SEA=reshape(SEAterms,size(Ahat,1),size(Ahat,2))';
                end
            end
            SEC=reshape(SECterms,size(Chat,2),size(Chat,1))';
            SEAlpha=reshape(SEAlphaterms,size(alphahat,1),size(alphahat,2));

            if(KFEM_Constraints.RhatDiag==1)
                SER=diag(SERterms);
            else
                SER=reshape(SERterms,size(Rhat,1),size(Rhat,2))';
            end
            if(KFEM_Constraints.QhatDiag==1)
                SEQ=diag(SEQterms);
            else
                SEQ=reshape(SEQterms,size(Qhat,1),size(Qhat,2))'; 
            end
            if(KFEM_Constraints.EstimateA==1)
                SE.A = SEA;
            end
            SE.Q = SEQ;
            SE.C = SEC;
            SE.R = SER;
            SE.alpha = SEAlpha;

            if(KFEM_Constraints.EstimatePx0==1)
                SE.Px0=SES;
            end
            if(KFEM_Constraints.Estimatex0==1)
                SE.x0=SEx0;
            end
            % Compute parameter p-values
            if(KFEM_Constraints.EstimateA==1)
                clear h p;
                if(KFEM_Constraints.AhatDiag==1)
                    VecParams = diag(Ahat);
                    VecSE     = diag(SEA);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pA = diag(p);
                else
                    VecParams = reshape(Ahat,[numel(Ahat) 1]);
                    VecSE     = reshape(SEA, [numel(Ahat) 1]);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end  
                    pA = reshape(p, [size(Ahat,1) size(Ahat,2)]);
                end
            end
            %C matrix
            clear h p;
            VecParams = reshape(Chat,[numel(Chat) 1]);
            VecSE     = reshape(SEC, [numel(Chat) 1]);
            for i=1:length(VecParams)
               [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end 
            pC = reshape(p, [size(Chat,1) size(Chat,2)]);

            %R matrix
            clear h p;
            if(KFEM_Constraints.RhatDiag==1)
                if(KFEM_Constraints.RhatIsotropic==1)
                    VecParams = Rhat(1,1);
                    VecSE     = SER(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pR = diag(p);
                else
                    VecParams = diag(Rhat);
                    VecSE     = diag(SER);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pR = diag(p);
                end
            else
                VecParams = reshape(Rhat,[numel(Rhat) 1]);
                VecSE     = reshape(SER, [numel(Rhat) 1]);
                for i=1:length(VecParams)
                   [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end  
                pR = reshape(p, [size(Rhat,1) size(Rhat,2)]);
            end

            %Q matrix
            clear h p;
            if(KFEM_Constraints.QhatDiag==1)
                if(KFEM_Constraints.QhatIsotropic==1)
                    VecParams = Qhat(1,1);
                    VecSE     = SEQ(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pQ = diag(p);
                else
                    VecParams = diag(Qhat);
                    VecSE     = diag(SEQ);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pQ = diag(p);
                end
            else
                VecParams = reshape(Qhat,[numel(Qhat) 1]);
                VecSE     = reshape(SEQ, [numel(Qhat) 1]);
                for i=1:length(VecParams)
                   [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end  
                pQ = reshape(p, [size(Qhat,1) size(Qhat,2)]);
            end
            %Px0
            if(KFEM_Constraints.EstimatePx0==1)
                clear h p;
                if(KFEM_Constraints.Px0Isotropic==1)
                    VecParams = Px0hat(1,1);
                    VecSE     = SES(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pPX0 = diag(p);
                else
                    VecParams = diag(Px0hat);
                    VecSE     = diag(SES);
                    for i=1:length(VecParams)
                        [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pPX0 = diag(p);
                end
            end

            clear h p;
            VecParams = alphahat;
            VecSE     = SEAlpha;
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end
            pAlpha = p';

            if(KFEM_Constraints.Estimatex0==1)
                clear h p;
                VecParams = x0hat;
                VecSE     = SEx0;
                for i=1:length(VecParams)
                    [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end
                pX0 = p';
            end
            if(KFEM_Constraints.EstimateA==1)
                Pvals.A = pA;
            end
            Pvals.Q = pQ;
            Pvals.C = pC;
            Pvals.R = pR;
            Pvals.alpha = pAlpha;
            if(KFEM_Constraints.EstimatePx0==1)
                Pvals.Px0 = pPX0;
            end
            if(KFEM_Constraints.Estimatex0==1)
                Pvals.x0 = pX0;
            end

        end
        function [x_K,W_K,logll,ExpectationSums]=KF_EStep(A,Q,C,R, y, alpha, x0, Px0)
             DEBUG = 0;

            Dx = size(A,2);
            Dy = size(C,1);
            K=size(y,2);
            [x_p, W_p, x_u, W_u] = DecodingAlgorithms.kalman_filter(A, C, Q, R,Px0, x0, y-alpha*ones(1,size(y,2)));
            
            [x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u);

            %Best estimates of initial states given the data
            W1G0 = A*Px0*A' + Q;
            L0=Px0*A'/W1G0;
            
            Ex0Gy = x0+L0*(x_K(:,1)-x_p(:,1));        
            Px0Gy = Px0+L0*(eye(size(W_K(:,:,1)))/(W_K(:,:,1))-eye(size(W1G0))/W1G0)*L0';
            Px0Gy = (Px0Gy+Px0Gy')/2;
            numStates = size(x_K,1);
            Wku=zeros(numStates,numStates,K,K);
            Tk = zeros(numStates,numStates,K-1);
            for k=1:K
                Wku(:,:,k,k)=W_K(:,:,k);
            end

            for u=K:-1:2
                for k=(u-1):-1:(u-1)
                    Tk(:,:,k)=A;
%                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
                    Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
                    Wku(:,:,u,k)=Wku(:,:,k,u)';
                end
            end
            
            %All terms
            Sxkm1xk = zeros(Dx,Dx);
            Sxkxkm1 = zeros(Dx,Dx);
            Sxkm1xkm1 = zeros(Dx,Dx);
            Sxkxk = zeros(Dx,Dx);
            Sykyk = zeros(Dy,Dy);
            Sxkyk = zeros(Dx,Dy);
            for k=1:K
                if(k==1)
                    Sxkm1xk   = Sxkm1xk+Px0*A'/W_p(:,:,1)*Wku(:,:,1,1);
                    Sxkm1xkm1 = Sxkm1xkm1+Px0+x0*x0';     
                else
                    Sxkm1xk =  Sxkm1xk+Wku(:,:,k-1,k)+x_K(:,k-1)*x_K(:,k)';
                    Sxkm1xkm1= Sxkm1xkm1+Wku(:,:,k-1,k-1)+x_K(:,k-1)*x_K(:,k-1)';
                end
                Sxkxk = Sxkxk+Wku(:,:,k,k)+x_K(:,k)*x_K(:,k)';
                Sykyk = Sykyk+(y(:,k)-alpha)*(y(:,k)-alpha)';
                Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alpha)';

            end
            Sx0x0 = Px0+x0*x0';
            Sxkxk = 0.5*(Sxkxk+Sxkxk');
            Sykyk = 0.5*(Sykyk+Sykyk');
            sumXkTerms = Sxkxk-A*Sxkm1xk-Sxkm1xk'*A'+A*Sxkm1xkm1*A';
            sumYkTerms = Sykyk - C*Sxkyk - Sxkyk'*C' + C*Sxkxk*C';      
            Sxkxkm1 = Sxkm1xk';
            
           

            logll = -Dx*K/2*log(2*pi)-K/2*log(det(Q))-Dy*K/2*log(2*pi)...
                    -K/2*log(det(R))- Dx/2*log(2*pi) -1/2*log(det(Px0))  ...
                    -1/2*trace((eye(size(Q))/Q)*sumXkTerms) ...
                    -1/2*trace((eye(size(R))/R)*sumYkTerms) ...
                    -Dx/2;
                string0 = ['logll: ' num2str(logll)];
                disp(string0);
                if(DEBUG==1)
                    string1 = ['-K/2*log(det(Q)):' num2str(-K/2*log(det(Q)))];
                    string11 = ['-K/2*log(det(R)):' num2str(-K/2*log(det(R)))];
                    string12= ['Constants: ' num2str(-Dx*K/2*log(2*pi)-Dy*K/2*log(2*pi)- Dx/2*log(2*pi) -Dx/2 -1/2*log(det(Px0)))];
                    string3 = ['-.5*trace(Q\sumXkTerms): ' num2str(-.5*trace(Q\sumXkTerms))];
                    string4 = ['-.5*trace(R\sumYkTerms): ' num2str(-.5*trace(R\sumYkTerms))];

                    disp(string1);
                    disp(['Q=' num2str(diag(Q)')]);
                    disp(string11);
                    disp(['R=' num2str(diag(R)')]);
                    disp(string12);
                    disp(string3);
                    disp(string4);
                end

                ExpectationSums.Sxkm1xkm1=Sxkm1xkm1;
                ExpectationSums.Sxkm1xk=Sxkm1xk;
                ExpectationSums.Sxkxkm1=Sxkxkm1;
                ExpectationSums.Sxkxk=Sxkxk;
                ExpectationSums.Sxkyk=Sxkyk;
                ExpectationSums.Sykyk=Sykyk;
                ExpectationSums.sumXkTerms=sumXkTerms;
                ExpectationSums.sumYkTerms=sumYkTerms;
                ExpectationSums.Sx0 = Ex0Gy;
                ExpectationSums.Sx0x0 = Px0Gy + Ex0Gy*Ex0Gy';

        end
        function [Ahat, Qhat, Chat, Rhat, alphahat, x0hat, Px0hat] = KF_MStep(y,x_K,x0, Px0, ExpectationSums,KFEM_Constraints)
            if(nargin<6 || isempty(KFEM_Constraints))
                KFEM_Constraints = DecodingAlgorithms.KF_EMCreateConstraints;
            end
            Sxkm1xkm1=ExpectationSums.Sxkm1xkm1;
            Sxkxkm1=ExpectationSums.Sxkxkm1;
            Sxkxk=ExpectationSums.Sxkxk;
            Sxkyk=ExpectationSums.Sxkyk;
            sumXkTerms = ExpectationSums.sumXkTerms;
            sumYkTerms = ExpectationSums.sumYkTerms;
            Sx0 = ExpectationSums.Sx0;
            Sx0x0 = ExpectationSums.Sx0x0;

            [N,K] = size(x_K); 
            
            
            if(KFEM_Constraints.AhatDiag==1)
                I=eye(N,N);
                Ahat = (Sxkxkm1.*I)/(Sxkm1xkm1.*I);
            else
                Ahat = Sxkxkm1/Sxkm1xkm1;
            end
            
            
%              [V,D] = eig(Px0hat);
%              D(D<0)=1e-9;
%              Px0hat = V*D*V';
            
            
             Chat = Sxkyk'/Sxkxk;             
             alphahat = sum(y - Chat*x_K,2)/K;
             
             if(KFEM_Constraints.QhatDiag==1)
                 if(KFEM_Constraints.QhatIsotropic==1)
                     Qhat=1/(N*K)*trace(sumXkTerms)*eye(N,N);
                 else
                     I=eye(N,N);
                     Qhat=1/K*(sumXkTerms.*I);
                     Qhat = (Qhat + Qhat')/2;
                 end
             else
                 Qhat=1/K*sumXkTerms;
                 Qhat = (Qhat + Qhat')/2;
             end
             dy=size(sumYkTerms,1);
             if(KFEM_Constraints.RhatDiag==1)
                 if(KFEM_Constraints.RhatIsotropic==1)
                     I=eye(dy,dy);
                     Rhat = 1/(dy*K)*trace(sumYkTerms)*I;
                 else
                     
                     I=eye(dy,dy);
                     Rhat = 1/K*(sumYkTerms.*I);
                     Rhat = (Rhat + Rhat')/2;
                 end
             else
                 Rhat = 1/K*(sumYkTerms);
                 Rhat = (Rhat + Rhat')/2;  
             end
             if(KFEM_Constraints.Estimatex0)
                x0hat = (inv(Px0)+Ahat'/Qhat*Ahat)\(Ahat'/Qhat*x_K(:,1)+Px0\x0);
            else
                x0hat = x0;
            end
             
            if(KFEM_Constraints.EstimatePx0==1)
                if(KFEM_Constraints.Px0Isotropic==1)
                   Px0hat=(trace(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0'))/(N*K))*eye(N,N); 
                else
                    I=eye(N,N);
                    Px0hat =(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0')).*I;
                    Px0hat = (Px0hat+Px0hat')/2;
                    [V,Lambda]=eig(Px0hat);
                    Lambda = diag(Lambda);
                    if(min(Lambda)<eps)
                        Lambda(Lambda==min(Lambda))=eps;
                        Px0hat = V*diag(Lambda)*V';
                    end
                    
                end
                
            else
                Px0hat =Px0;
            end
            
        end
        
        %% Mixed Point Process and Continuous Observation (mPPCO)
        function [x_pLag, W_pLag, x_uLag, W_uLag] = mPPCO_fixedIntervalSmoother(A, Q, C, R, y, alpha, dN,lags,mu,beta,fitType,delta,gamma,windowTimes,x0,Px0,HkAll)    
            nStates = size(A,2);

            [numCells,N]   = size(dN); % N time samples
            nObs = size(C,1);
            ns=size(A,1); % number of states

            if(nargin<17 || isempty(HkAll))
                HkAll=[];
            end
            if(nargin<16 || isempty(Px0))
                Px0 = zeros(ns,ns);
            end
            
            if(nargin<15 || isempty(x0))
               x0=zeros(ns,1);

            end
            if(nargin<14 || isempty(windowTimes))
               windowTimes=[]; 
            end
            if(nargin<13 || isempty(gamma))
                gamma=0;
            end
            if(nargin<12 || isempty(delta))
                delta = .001;
            end

            
            minTime=0;
            maxTime=(size(dN,2)-1)*delta;

            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                HkAll = zeros(size(dN,2),length(windowTimes)-1,numCells);
                for c=1:numCells
                    nst{c} = nspikeTrain( (find(dN(c,:)==1)-1)*delta);
                    nst{c}.setMinTime(minTime);
                    nst{c}.setMaxTime(maxTime);
                    nst{c}=nst{c}.resample(1/delta);
                    HkAll(:,:,c) = histObj.computeHistory(nst{c}).dataToMatrix;
    %                 HkAll{c} = histObj.computeHistory(nst{c}).dataToMatrix;
                end
                if(size(gamma,2)==1 && numCells>1) % if more than 1 cell but only 1 gamma
                    gammaNew(:,c) = gamma;
                else
                    gammaNew=gamma;
                end
                gamma = gammaNew;

            else
                for c=1:numCells
    %                 HkAll{c} = zeros(N,1);
                    HkAll(:,:,c) = zeros(N,1);
                    gammaNew(c)=0;
                end
                gamma=gammaNew;

            end
            if(size(gamma,2)~=numCells)
                gamma=gamma';
            end
        
                
            Alag = zeros((lags+1)*nStates,(lags+1)*nStates,N);
            Qlag = zeros((lags+1)*nStates,(lags+1)*nStates,N);
            Clag = zeros(nObs,(lags+1)*nStates,N);
            Rlag = zeros(nObs,nObs,N);
            x0lag = zeros(length(x0)*(lags+1),1);
            Px0lag = zeros((lags+1)*nStates,(lags+1)*nStates);
            Px0lag((1:nStates),(1:nStates))=Px0;
            x0lag(1:nStates,1)=x0;
            for n=1:N
                offset = 0;
                for i=1:(lags+1)
                    if(i==1)
                        Alag((1:nStates)+offset,(1:nStates)+offset,n)=A(:,:,min(size(A,3),n));
                        Qlag((1:nStates)+offset,(1:nStates)+offset,n)=Q(:,:,min(size(Q,3),n));
                        Clag((1:nObs),(1:nStates)+offset,n)=C(:,:,min(size(C,3),n));
                        Rlag((1:nObs),(1:nObs),n) = R(:,:,min(size(R,3),n));

                    else
                        Alag((1:nStates)+offset,(1:nStates)+(offset-nStates),n)=eye(nStates,nStates);
                        Qlag((1:nStates)+offset,(1:nStates)+offset,n)=zeros(nStates,nStates);
                        Clag((1:nObs),(1:nStates)+offset,n)=zeros(nObs,nStates);

                    end
                    offset=offset+nStates;
                end
            end
            
            betaLag = zeros((lags+1)*nStates, numCells);
            betaLag(1:nStates,1:numCells)=beta;
            [x_p, W_p, x_u, W_u] = DecodingAlgorithms.mPPCODecodeLinear(Alag, Qlag, Clag, Rlag, y, alpha, dN,mu,betaLag,fitType,delta,gamma,windowTimes,x0lag,Px0lag,HkAll);
            

            x_pLag = x_p((lags*nStates+1):(lags+1)*nStates,:);
            W_pLag = W_p((lags*nStates+1):(lags+1)*nStates,(lags*nStates+1):(lags+1)*nStates,:);
            x_uLag = x_u((lags*nStates+1):(lags+1)*nStates,:);
            W_uLag = W_u((lags*nStates+1):(lags+1)*nStates,(lags*nStates+1):(lags+1)*nStates,:);
           
        end
        function [x_p, W_p, x_u, W_u] = mPPCODecodeLinear(A, Q, C, R, y, alpha, dN,mu,beta,fitType,delta,gamma,windowTimes,x0,Px0,HkAll)
        % [x_p, W_p, x_u, W_u] = mPPCODecodeLinear(A, Q, C, R, y, dN, mu, beta,fitType, delta, gamma,windowTimes, x0)
        % Point process adaptive filter with the assumption of linear
        % expresion for the conditional intensity functions (see below). If
        % the terms in the conditional intensity function include
        % polynomial powers of a variable for example, these expressions do
        % not hold. Use the PPDecodeFilter instead since it will compute
        % these expressions symbolically. However, because of the matlab
        % symbolic toolbox, it runs much slower than this version.
        % 
        % Assumes in both cases that 
        %   x_t = A*x_{t-1} + v_{t}     w_{t} ~ Normal with zero me and
        %                                       covariance Q
        %
        %   y_t = C*x_{t} + w_{t}     w_{t} ~ Normal with zero me and
        %                                       covariance R
        %
        % Paramerters:
        %  
        % A:        The state transition matrix from the x_{t-1} to x_{t}
        %
        % Q:        The covariance of the process noise v_t
        %
        % C:        The observation matrix
        %
        % R:        The covariance of the observation noise w_t
        %
        % y:        The continuous observations
        %
        % alpha:    Offset for the observations
        %
        % dN:       A C x N matrix of ones and zeros corresponding to the
        %           observed spike trains. N is the number of time steps in
        %           my code. C is the number of cells
        %
        % mu:       Cx1 vector of baseline firing rates for each cell. In
        %           the CIF expression in 'fitType' description 
        %           mu_c=mu(c);
        %
        % beta:     nsxC matrix of coefficients for the conditional
        %           intensity function. ns is the number of states in x_t 
        %           In the conditional intesity function description below
        %           beta_c = beta(:,c)';
        %
        % fitType: 'poisson' or 'binomial'. Determines how the beta and
        %           gamma coefficients are used to compute the conditional
        %           intensity function.
        %           For the cth cell:
        %           If poisson: lambda*delta = exp(mu_c+beta_c*x + gamma_c*hist_c)
        %           If binomial: logit(lambda*delta) = mu_c+beta_c*x + gamma_c*hist_c
        %
        % delta:    The number of seconds per time step. This is used to compute
        %           th history effect for each spike train using the input
        %           windowTimes and gamma
        %
        % gamma:    length(windowTimes)-1 x C matrix of the history
        %           coefficients for each window in windowTimes. In the 'fitType'
        %           expression above:
        %           gamma_c = gamma(:,c)';
        %           If gamma is a length(windowTimes)-1x1 vector, then the
        %           same history coefficients are used for each cell.
        %
        % windowTimes: Defines the distinct windows of time (in seconds)
        %           that will be computed for each spike train.
        %
        % x0:       The initial state
        %
        % Px0:      The initial state covariance
        
        
        
        [numCells,N]   = size(dN); % N time samples, C cells
        ns=size(A,1); % number of states
        if(nargin<16 || isempty(HkAll))
            HkAll=[];
        end
        if(nargin<15 || isempty(Px0))
           Px0=zeros(ns,ns);
        end
        if(nargin<14 || isempty(x0))
           x0=zeros(ns,1);
           
        end
        if(nargin<13 || isempty(windowTimes))
           windowTimes=[]; 
        end
        if(nargin<12 || isempty(gamma))
            gamma=0;
        end
        if(nargin<11 || isempty(delta))
            delta = .001;
        end
        
        
        minTime=0;
        maxTime=(size(dN,2)-1)*delta;
        
%         numCells=size(dN,1);
        if(~isempty(HkAll))
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for c=1:numCells
                    nst{c} = nspikeTrain( (find(dN(c,:)==1)-1)*delta);
                    nst{c}.setMinTime(minTime);
                    nst{c}.setMaxTime(maxTime);
                    nst{c}=nst{c}.resample(1/delta);
%                     HkAll{c} = histObj.computeHistory(nst{c}).dataToMatrix;
                    HkAll(:,:,c) = histObj.computeHistory(nst{c}).dataToMatrix;
                end
                if(size(gamma,2)==1 && numCells>1) % if more than 1 cell but only 1 gamma
                    gammaNew(:,c) = gamma;
                else
                    gammaNew = gamma;
                end
                gamma = gammaNew;
            end

        else
            for c=1:numCells
%                 HkAll{c} = zeros(N,1);
                HkAll(:,:,c) = zeros(N,1);
                gammaNew(c)=0;
            end
            gamma=gammaNew;

        end
        

        
        %% Initialize the numCells
        x_p     = zeros( size(A,2), N+1 );
        x_u     = zeros( size(A,2), N );
        W_p    = zeros( size(A,2),size(A,2), N+1 );
        W_u    = zeros( size(A,2),size(A,2), N );
        A1=A(:,:,min(size(A,3),1));
        x_p(:,1) = A1*x0;
        W_p(:,:,1) = A1 * Px0 * A1' +Q(:,:,min(size(Q,3),1));
        Histtermperm = permute(HkAll,[2 3 1]);
%         WuConv = [];
        for n=1:N
%             [x_u, W_u,lambdaDeltaMat] = mPPCODecode_update(x_p, W_p, C, R, y, alpha, dN,mu,beta,fitType,gamma,HkAll,time_index,WuConv)
            [x_u(:,n), W_u(:,:,n)] = DecodingAlgorithms.mPPCODecode_update(x_p(:,n), W_p(:,:,n),  C(:,:,min(size(C,3),n)), R(:,:,min(size(R,3),n)), y(:,n), alpha(:,min(size(alpha,3),n)),dN,mu,beta,fitType,gamma,Histtermperm,n,[]); %expects History with time on 3rd index
            if(n<N)
                [x_p(:,n+1), W_p(:,:,n+1)] = DecodingAlgorithms.mPPCODecode_predict(x_u(:,n), W_u(:,:,n), A(:,:,min(size(A,3),n)), Q(:,:,min(size(Q,3),n)));
            end
%             if(n>1 && isempty(WuConv))
%                 diffWu = abs(W_u(:,:,n)-W_u(:,:,n-1));
%                 maxWu  = max(max(diffWu));
%                 if(maxWu<5e-4)
%                     WuConv = W_u(:,:,n);
%                     WuConvIter = n;
%                 end
%             end
        end
      
        
        end
        function [x_p, W_p] = mPPCODecode_predict(x_u, W_u, A, Q)
            x_p     = A * x_u;
            W_p    = A * W_u * A' + Q;
%             if(rcond(W_p)<1000*eps)
%                 W_p=W_u; % See Srinivasan et al. 2007 pg. 529
%             end
            W_p = .5*(W_p + W_p'); %To help with symmetry of matrix;

        end 
        function [x_u, W_u,lambdaDeltaMat] = mPPCODecode_update(x_p, W_p, C, R, y, alpha, dN,mu,beta,fitType,gamma,HkAll,time_index,WuConv)
            [numCells,N]   = size(dN); % N time samples, C cells
            if(nargin<13 || isempty(WuConv))
                WuConv=[];
            end
            if(nargin<12 || isempty(time_index))
                time_index=1;
            end
            if(nargin<11 || isempty(HkAll))
                  HkAll = zeros(numCells,1);
            end
            if(nargin<10 || isempty(gamma))
                gamma=zeros(1,numCells);
            end
            if(nargin<9 || isempty(fitType))
                fitType = 'poisson';
            end


            sumValVec=zeros(size(W_p,1),1);
            sumValMat=zeros(size(W_p,2),size(W_p,2));
            lambdaDeltaMat = zeros(numCells,1);

            if(numel(gamma)==1 && gamma==0)
                gamma = zeros(size(mu))';
            end
            if(strcmp(fitType,'binomial'))
                Histterm = HkAll(:,:,time_index);
                if(size(Histterm,1)~=numCells) %make sure Histterm has proper orientation
                    Histterm = Histterm';
                end

                if(size(gamma,2)~=size(mu,1)) 
                    if(size(gamma,1)==size(Histterm,1)) %All cells have same history
                        gamma = repmat(gamma,[1 numCells]);
                    end
                end
                    linTerm = mu+beta'*x_p + diag(gamma'*Histterm');
                    lambdaDeltaMat = exp(linTerm)./(1+exp(linTerm));
                    if(any(isnan(lambdaDeltaMat))||any(isinf(lambdaDeltaMat)))
                        indNan = isnan(lambdaDeltaMat);
                        indInf = isinf(lambdaDeltaMat);
                        lambdaDeltaMat(indNan)=1;
                        lambdaDeltaMat(indInf)=1;
                    end
                    sumValVec=sum(repmat(((dN(:,time_index)-lambdaDeltaMat(:,1)).*(1-lambdaDeltaMat(:,1)))',size(beta,1),1).*beta,2);
                    tempVec = ((dN(:,time_index)+(1-2*(lambdaDeltaMat(:,1)))).*(1-(lambdaDeltaMat(:,1))).*(lambdaDeltaMat(:,1)))';
%                     tempVec((tempVec<0))=0;
%                     tempVec((tempVec>1))=1;
                    sumValMat = (repmat(tempVec,size(beta,1),1).*beta)*beta';
            elseif(strcmp(fitType,'poisson'))
                Histterm = HkAll(:,:,time_index);
                if(size(Histterm,1)~=numCells) %make sure Histterm has proper orientation
                    Histterm = Histterm';
                end

                if(size(gamma,2)~=size(mu,1)) 
                    if(size(gamma,1)==size(Histterm,1)) %All cells have same history
                        gamma = repmat(gamma,[1 numCells]);
                    end
                end
                        
                linTerm = mu+beta'*x_p + diag(gamma'*Histterm');
                lambdaDeltaMat = exp(linTerm);
                if(any(isnan(lambdaDeltaMat))||any(isinf(lambdaDeltaMat)))
                    indNan = isnan(lambdaDeltaMat);
                    indInf = isinf(lambdaDeltaMat);
                    lambdaDeltaMat(indNan)=1;
                    lambdaDeltaMat(indInf)=1;
                end
                sumValVec=sum(repmat(((dN(:,time_index)-lambdaDeltaMat(:,1)))',size(beta,1),1).*beta,2);
                sumValMat = (repmat(lambdaDeltaMat(:,1)',size(beta,1),1).*beta)*beta';
            end
            if(isempty(WuConv))
                sumValMat = sumValMat+C'*(R\C);
                I=eye(size(W_p));
                Wu=W_p*(I-(I+sumValMat*W_p)\(sumValMat*W_p));
                if(any(any(isnan(Wu)))||any(any(isinf(Wu))))
                    Wu=W_p;
                end
               % Make sure that the update covariance is positive definite.
                W_u = Wu;
                W_u = .5*(W_u + W_u'); %To help with symmetry of matrix;
            else
                W_u = WuConv;
            end
            x_u     = x_p + W_u*(sumValVec)+((W_u*C')/R)*(y-C*x_p -alpha);


        end
        function C = mPPCO_EMCreateConstraints(EstimateA, AhatDiag,QhatDiag,QhatIsotropic,RhatDiag,RhatIsotropic,Estimatex0,EstimatePx0, Px0Isotropic,mcIter,EnableIkeda)
            %By default, all parameters are estimated. To empose diagonal
            %structure on the EM parameter results must pass in the
            %constraints element
            if(nargin<11 || isempty(EnableIkeda))
                EnableIkeda=0;
            end
            if(nargin<10 || isempty(mcIter))
                mcIter=1000;
            end
            if(nargin<9 || isempty(Px0Isotropic))
                Px0Isotropic=0;
            end
            if(nargin<8 || isempty(EstimatePx0))
                EstimatePx0=1;
            end
            if(nargin<7 || isempty(Estimatex0))
                Estimatex0=1;
            end
            if(nargin<6 || isempty(RhatIsotropic))
                RhatIsotropic=0;
            end
            if(nargin<5 || isempty(RhatDiag))
                RhatDiag=1;
            end
            if(nargin<4 || isempty(QhatIsotropic))
                QhatIsotropic=0;
            end
            if(nargin<3 || isempty(QhatDiag))
                QhatDiag=1;
            end
            if(nargin<2)
                AhatDiag=0;
            end
            if(nargin<1)
                EstimateA=1;
            end
            C.EstimateA= EstimateA;
            C.AhatDiag = AhatDiag;
            C.QhatDiag = QhatDiag;
            if(QhatDiag && QhatIsotropic)
                C.QhatIsotropic=1;
            else
                C.QhatIsotropic=0;
            end
            C.RhatDiag = RhatDiag;
            if(RhatDiag && RhatIsotropic)
                C.RhatIsotropic=1;
            else
                C.RhatIsotropic=0;
            end
            C.Estimatex0 = Estimatex0;
            C.EstimatePx0 = EstimatePx0;
            if(EstimatePx0 && Px0Isotropic)
                C.Px0Isotropic=1;
            else
                C.Px0Isotropic=0; 
            end
            C.mcIter = mcIter;
            C.EnableIkeda = EnableIkeda;
        end  
        function [SE,Pvals,nTerms] = mPPCO_ComputeParamStandardErrors(y, dN, xKFinal, WKFinal, Ahat, Qhat, Chat, Rhat, alphahat, x0hat, Px0hat, ExpectationSumsFinal, fitType, muhat, betahat, gammahat, windowTimes, HkAll, mPPCOEM_Constraints)

            % Use inverse observed information matrix to estimate the standard errors of the estimated model parameters
         % Requires computation of the complete information matrix and an estimate of the missing information matrix

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Complete Information Matrices     
        % Recall from McLachlan and Krishnan Eq. 4.7
        %    Io(theta;y) = Ic(theta;y) - Im(theta;y)
        %    Io(theta;y) = Ic(theta;y) - cov(Sc(X;theta)Sc(X;theta)')
        % where Sc(X;theta) is the score vector of the complete log likelihood
        % function evaluated at theta. We first compute Ic term by term and then
        % approximate the covariance term using Monte Carlo approximation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(nargin<19 || isempty(mPPCOEM_Constraints))
                mPPCOEM_Constraints=DecodingAlgorithms.mPPCO_EMCreateConstraints;
            end

            if(mPPCOEM_Constraints.EstimateA==1)
                if(mPPCOEM_Constraints.AhatDiag==1)
                    IAComp=zeros(numel(diag(Ahat)),numel(diag(Ahat)));
                else
                    IAComp=zeros(numel(Ahat),numel(Ahat));
                end
                [n1,n2] =size(Ahat);
                el=(eye(n1,n1));
                em=(eye(n2,n2));
                cnt=1;
                N=size(y,2);

                if(mPPCOEM_Constraints.AhatDiag==1)
                    for l=1:n1
                        for m=l
                            termMat=Qhat\el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1.*eye(n1,n2);
                            termvec = diag(termMat);
                            IAComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                else
                    for l=1:n1
                        for m=1:n2
                            termMat=(inv(Qhat))*el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1;
                            termvec=reshape(termMat',1,numel(Ahat));
                            IAComp(:,cnt)=termvec';
                            cnt=cnt+1;
                        end
                    end
                end
            end


            ICComp=zeros(numel(Chat),numel(Chat));
            [n1,n2] =size(Chat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            for l=1:n1
                for m=1:n2
                    termMat=Rhat\el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkxk;
                    termvec=reshape(termMat',1,numel(Chat));
                    ICComp(:,cnt)=termvec';
                    cnt=cnt+1;
                end
            end

            [n1,n2] =size(Rhat);
            ei=(eye(n1,n1));
            ej=(eye(n2,n2));
            cnt=1;
            [dy,N]=size(y);
            dx=size(xKFinal,1);

            [n1,n2] =size(Rhat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            N=size(y,2);
            if(mPPCOEM_Constraints.RhatDiag==1)
                if(mPPCOEM_Constraints.RhatIsotropic==1)
                    IRComp = 0.5*N*dy*Rhat(1,1)^(-2);
                else
                    IRComp=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
                    for l=1:n1
                        for m=l
                            termMat= N/2*(Rhat)\em(:,m)*el(:,l)'/(Rhat);
                            termvec=diag(termMat);
                            IRComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            else
                IRComp=zeros(numel(diag(Rhat)),numel(diag(Rhat)));
                for l=1:n1
                    for m=1:n2
                        termMat= N/2*(Rhat)\em(:,m)*el(:,l)'/(Rhat);
                        termvec=reshape(termMat',1,numel(Rhat));
                        IRComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            end

            [n1,n2] =size(Qhat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            N=size(y,2);
            if(mPPCOEM_Constraints.QhatDiag==1)
                if(mPPCOEM_Constraints.QhatIsotropic==1)
                    IQComp=zeros(1,1);
                    IQComp =  0.5*N*dx*Qhat(1,1)^(-2); 
                else
                    IQComp=zeros(numel(diag(Qhat)),numel(diag(Qhat)));
                    for l=1:n1
                        for m=l
                            termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                            termvec=diag(termMat);
                            IQComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            else
                IQComp=zeros(numel(Qhat),numel(Qhat));
                for l=1:n1
                    for m=1:n2
                        termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                        termvec=reshape(termMat',1,numel(Qhat));
                        IQComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            end

            if(mPPCOEM_Constraints.EstimatePx0==1)
                if(mPPCOEM_Constraints.Px0Isotropic==1)
                    ISComp =  0.5*dx*Px0hat(1,1)^(-2);
                else
                    ISComp=zeros(numel(diag(Px0hat)),numel(diag(Px0hat)));
                    [n1,n2] =size(Px0hat);
                    el=(eye(n1,n1));
                    em=(eye(n2,n2));
                    cnt=1;
                    for l=1:n1
                        for m=l
                            termMat= 1/2*(Px0hat)\em(:,m)*el(:,l)'/(Px0hat);
                            termvec=diag(termMat);
                            ISComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            end

            if(mPPCOEM_Constraints.Estimatex0==1)
                Ix0Comp=eye(size(Px0hat))/Px0hat+(Ahat'/Qhat)*Ahat;
            end

            IAlphaComp = N*eye(size(Rhat))/Rhat;
            K=size(y,2);
            numCells=size(betahat,2);
%             McExp=500;    
            McExp=mPPCOEM_Constraints.mcIter;
            xKDrawExp = zeros(size(xKFinal,1),K,McExp);
            

            % Generate the Monte Carlo
            for k=1:K
%                 WuTemp=squeeze(WKFinal(:,:,k));
                WuTemp=(WKFinal(:,:,k));
                [chol_m,p]=chol(WuTemp);
                z=normrnd(0,1,size(xKFinal,1),McExp);
                xKDrawExp(:,k,:)=repmat(xKFinal(:,k),[1 McExp])+(chol_m*z);
            end
            
            IBetaComp =zeros(size(xKFinal,1)*numCells,size(xKFinal,1)*numCells);
            xkPerm = permute(xKDrawExp,[1 3 2]);
            pools = matlabpool('size'); %number of parallel workers
            if(pools==0)
                if(strcmp(fitType,'poisson'))
                    for c=1:numCells
                        HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1));
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            Wk = WKFinal(:,:,k);

    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk=xkPerm(:,:,k);
                           if(size(Hk,1)==numCells)
                               Hk = Hk';
                           end

                            if(numel(gammahat)==1)
                                gammaC=gammahat;
    %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat(:,c);
                            end

                            terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld=exp(terms);

                            HessianTerm=HessianTerm-1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                        end
                        startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                        IBetaComp(startInd:endInd,startInd:endInd)=-HessianTerm;
                    end
                else
                    for c=1:numCells
                        HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1));
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            Wk = WKFinal(:,:,k);
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = (xkPerm(:,:,k));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end

                            if(numel(gammahat)==1)
                                gammaC=gammahat;
    %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld=exp(terms)./(1+exp(terms));
                            ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                            ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                            ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                            HessianTerm=HessianTerm+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;

                        end
                        startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                        IBetaComp(startInd:endInd,startInd:endInd)=-HessianTerm;
                    end
                end
            else
                if(strcmp(fitType,'poisson'))
                    for c=1:numCells
                        HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1),K);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            Wk = WKFinal(:,:,k);

    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk=xkPerm(:,:,k);
                           if(size(Hk,1)==numCells)
                               Hk = Hk';
                           end

                            if(numel(gammahat)==1)
                                gammaC=gammahat;
    %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat(:,c);
                            end

                            terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld=exp(terms);

                            HessianTerm(:,:,k)=-1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                        end
                        startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                        IBetaComp(startInd:endInd,startInd:endInd)=-sum(HessianTerm,3);
                    end
                else
                    for c=1:numCells
                        HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1),K);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            Wk = WKFinal(:,:,k);
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = (xkPerm(:,:,k));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end

                            if(numel(gammahat)==1)
                                gammaC=gammahat;
    %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld=exp(terms)./(1+exp(terms));
                            ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                            ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                            ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                            HessianTerm(:,:,k)=+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;

                        end
                        startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                        IBetaComp(startInd:endInd,startInd:endInd)=-sum(HessianTerm,3);
                    end
                            
                end
            end
                        

            %CIF means
            IMuComp=zeros(numel(muhat),numel(muhat));
            xkPerm = permute(xKDrawExp,[1 3 2]);
            if(pools==0)
                for c=1:numCells
                    if(strcmp(fitType,'poisson'))
                        HessianTerm = 0;
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(:,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                            ld = exp(terms);
                            HessianTerm=HessianTerm-1/McExp*sum(ld,2);
                        end
                    elseif(strcmp(fitType,'binomial'))
                        HessianTerm = 0;
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(:,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                            ld = exp(terms)./(1+exp(terms));
                            ExplambdaDelta = 1/McExp*sum(ld,2);
                            ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                            ExplambdaDeltaCubed = 1/McExp*sum(ld.^3,2);
                            HessianTerm = HessianTerm -(dN(c,k)+1)*ExplambdaDelta ...
                                +(dN(c,k)+3)*ExplambdaDeltaSquare-3*ExplambdaDeltaCubed;
                        end
                    end
                    IMuComp(c,c) = -HessianTerm;
                end
            else
                for c=1:numCells
                    if(strcmp(fitType,'poisson'))
                        HessianTerm = zeros(K,1);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld = exp(terms);
                            HessianTerm(k)=-1/McExp*sum(ld,2);
                        end
                    elseif(strcmp(fitType,'binomial'))
                        HessianTerm = zeros(K,1);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ExplambdaDelta = 1/McExp*sum(ld,2);
                            ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                            ExplambdaDeltaCubed = 1/McExp*sum(ld.^3,2);
                            HessianTerm(k) =  -(dN(c,k)+1)*ExplambdaDelta ...
                                +(dN(c,k)+3)*ExplambdaDeltaSquare-3*ExplambdaDeltaCubed;
                        end
                    end
                    IMuComp(c,c) = -sum(HessianTerm);
                end
            end
            
            
            % Gamma Information Matrix
            IGammaComp = zeros(numel(gammahat),numel(gammahat));
            if(~isempty(windowTimes) && any(any(gammahat~=0)))
                xkPerm = permute(xKDrawExp,[1 3 2]);
                if(pools==0)
                     for c=1:numCells
                       if(strcmp(fitType,'poisson'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2));
                            for k=1:K
    %                             Hk = squeeze(HkAll(:,:,c));
                                Hk = (HkAll(:,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                ld = exp(terms);
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                HessianTerm=HessianTerm-Hk(k,:)'*Hk(k,:)*ExplambdaDelta;
                            end
                       elseif(strcmp(fitType,'binomial'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2));
                            for k=1:K
                                Hk = (HkAll(:,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                ld = exp(terms)./(1+exp(terms));
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                                ExplambdaDeltaCubed  = 1/McExp*sum(ld.^2,2);
                                HessianTerm=HessianTerm+(-ExplambdaDelta*(dN(c,k)+1)...
                                    +ExplambdaDeltaSquare*(dN(c,k)+3)...
                                    -2*ExplambdaDeltaCubed)*Hk(k,:)'*Hk(:,k);
                            end
                       end
                       startInd=size(HkAll,2)*(c-1)+1; endInd = size(HkAll,2)*c;
                       IGammaComp(startInd:endInd,startInd:endInd) = -HessianTerm;
                     end

                else
            
                    for c=1:numCells
                       if(strcmp(fitType,'poisson'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2),K);
                            parfor k=1:K
    %                             Hk = squeeze(HkAll(:,:,c));
                                Hk = (HkAll(k,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                                ld = exp(terms);
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                HessianTerm(:,:,k)=-Hk'*Hk*ExplambdaDelta;
                            end
                        elseif(strcmp(fitType,'binomial'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2),K);

                            parfor k=1:K
                                Hk = (HkAll(k,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                                ld = exp(terms)./(1+exp(terms));
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                                ExplambdaDeltaCubed  = 1/McExp*sum(ld.^2,2);
                                HessianTerm(:,:,k)=+(-ExplambdaDelta*(dN(c,k)+1)...
                                    +ExplambdaDeltaSquare*(dN(c,k)+3)...
                                    -2*ExplambdaDeltaCubed)*Hk'*Hk;
                            end
                       end
                       startInd=size(HkAll,2)*(c-1)+1; endInd = size(HkAll,2)*c;
                       IGammaComp(startInd:endInd,startInd:endInd) = -sum(HessianTerm,3);
                    end

                end
            end
        
              
            
            if(mPPCOEM_Constraints.EstimateA==1)
                n1=size(IAComp,1); 
            else
                n1=0;
            end
            n2=size(IQComp,1); n3=size(ICComp,1); n4=size(IRComp,1); 
            if(mPPCOEM_Constraints.EstimatePx0==1)
                n5=size(ISComp,1); 
            else
                n5=0;
            end
            if(mPPCOEM_Constraints.Estimatex0==1)   
                n6=size(Ix0Comp,1);
            else
                n6=0;
            end
            n7=size(IAlphaComp,1);
            n8=size(IMuComp,1);
            n9=size(IBetaComp,1);
            if(numel(gammahat)==1)
                if(gammahat==0)
                    n10=0;
                end
            else
                n10=size(IGammaComp,1);
            end
            nTerms=n1+n2+n3+n4+n5+n6+n7+n8+n9+n10;
            IComp = zeros(nTerms,nTerms);
            if(mPPCOEM_Constraints.EstimateA==1)
                IComp(1:n1,1:n1)=IAComp;
            end
            offset=n1+1;
            IComp(offset:(n1+n2),offset:(n1+n2))=IQComp;
            offset=n1+n2+1;
            IComp(offset:(n1+n2+n3),offset:(n1+n2+n3))=ICComp;
            offset=n1+n2+n3+1;
            IComp(offset:(n1+n2+n3+n4),offset:(n1+n2+n3+n4))=IRComp;
            offset=n1+n2+n3+n4+1;
            if(mPPCOEM_Constraints.EstimatePx0==1);
                IComp(offset:(n1+n2+n3+n4+n5),offset:(n1+n2+n3+n4+n5))=ISComp;
            end
            offset=n1+n2+n3+n4+n5+1;
            if(mPPCOEM_Constraints.Estimatex0==1)
                IComp(offset:(n1+n2+n3+n4+n5+n6),offset:(n1+n2+n3+n4+n5+n6))=Ix0Comp;
            end
            offset=n1+n2+n3+n4+n5+n6+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7),offset:(n1+n2+n3+n4+n5+n6+n7))=IAlphaComp;
            offset=n1+n2+n3+n4+n5+n6+n7+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7+n8),offset:(n1+n2+n3+n4+n5+n6+n7+n8))=IMuComp;
            offset=n1+n2+n3+n4+n5+n6+n7+n8+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7+n8+n9),offset:(n1+n2+n3+n4+n5+n6+n7+n8+n9))=IBetaComp;
            offset=n1+n2+n3+n4+n5+n6+n7+n8+n9+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10),offset:(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10))=IGammaComp;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Missing Information Matrix
            %Approximate cov(Sc(X;theta)Sc(X;theta)')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Mc=mPPCOEM_Constraints.mcIter;
            xKDraw = zeros(size(xKFinal,1),N,Mc);

            % Generate the Monte Carlo samples for the unobserved data
            for n=1:N
                WuTemp=(WKFinal(:,:,n));
                [chol_m,p]=chol(WuTemp);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                xKDraw(:,n,:)=repmat(xKFinal(:,n),[1 Mc])+(chol_m*z);
            end


            if(mPPCOEM_Constraints.EstimatePx0|| mPPCOEM_Constraints.Estimatex0)
                [chol_m,p]=chol(Px0hat);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                x0Draw=repmat(x0hat,[1 Mc])+(chol_m*z); 
            else
               x0Draw=repmat(x0hat, [1 Mc]);

            end

            IMc = zeros(nTerms,nTerms,Mc);
            % Emperically estimate the covariance of the score
            pools = matlabpool('size'); %number of parallel workers 
            if(pools==0) % parallel toolbox is not enabled;
                for c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Dy=size(y,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);
                    Sykyk = zeros(Dy,Dy);
                    Sxkyk = zeros(Dx,Dy);        

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                        Sykyk = Sykyk+(y(:,k)-alphahat)*(y(:,k)-alphahat)';
                        Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alphahat)';

                    end
                    Sx0x0 = x_0*x_0';
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    Sykyk = 0.5*(Sykyk+Sykyk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    sumYkTerms = Sykyk - Chat*Sxkyk - Sxkyk'*Chat' + Chat*Sxkxk*Chat';      
                    Sxkxkm1 = Sxkm1xk';
                    Sykxk = Sxkyk';

                    

                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    sumYkTerms=0.5*(sumYkTerms+sumYkTerms');
                    if(mPPCOEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(mPPCOEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    else
                        ScoreAMc=[];
                    end

                    ScorC=Rhat\(Sykxk-Chat*Sxkxk);
                    ScoreCMc=reshape(ScorC',numel(ScorC),1);

                    if(mPPCOEM_Constraints.QhatDiag)
                        if(mPPCOEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end


                    ScoreAlphaMc = sum(Rhat\(y-Chat*x_K-alphahat*ones(1,N)),2);
                    if(mPPCOEM_Constraints.RhatDiag)
                        if(mPPCOEM_Constraints.RhatIsotropic)
                            ScoreR  =-.5*(K*Dy*Rhat(1,1)^(-1) - Rhat(1,1)^(-2)*trace(sumYkTerms));
                        else
                            ScoreR  =(-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat)));
                        end
                        ScoreRMc = diag(ScoreR);
                    else
                        ScoreR   =-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat));
                        ScoreRMc =reshape(ScoreR',numel(ScoreR),1);
                    end


                    if(mPPCOEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);
                    ScoreMuMc=zeros(numCells,1);
                    ScoreBetaMc=[];
                    ScoreGammaMc=[];
                    % Cell Scores
                    for nc=1:numCells
                        if(strcmp(fitType,'poisson'))
                            Hk = (HkAll(:,:,nc));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms);
                            ScoreMuMc(nc) = sum(dN(nc,:)-ld,2);
                            ScoreBetaMc = [ScoreBetaMc; sum(repmat((dN(nc,:)-ld),[Dx 1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-ld,[nHist 1]).*Hk',2)];
                        elseif(strcmp(fitType,'binomial'))
                            Hk = (HkAll(:,:,nc));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
                            nHist = size(Hk,2);
                            
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ScoreMuMc(nc) = sum(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,2);
                            ScoreBetaMc = [ScoreBetaMc;sum(repmat(dN(nc,:).*(1-ld) - ld.*(1-ld),[Dx,1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,[nHist 1]).*Hk',2)];
                        end
                        
                    end
                    ScoreVec = [ScoreAMc; ScoreQMc; ScoreCMc; ScoreRMc];
                    if(mPPCOEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(mPPCOEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreAlphaMc];
                    ScoreVec = [ScoreVec; ScoreMuMc; ScoreBetaMc];
                    if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                        ScoreVec=[ScoreVec;ScoreGammaMc];
                    end
                    
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end
            else %Use the parallel toolbox
                parfor c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Dy=size(y,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);
                    Sykyk = zeros(Dy,Dy);
                    Sxkyk = zeros(Dx,Dy);        

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                        Sykyk = Sykyk+(y(:,k)-alphahat)*(y(:,k)-alphahat)';
                        Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alphahat)';

                    end
                    Sx0x0 = x_0*x_0';
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    Sykyk = 0.5*(Sykyk+Sykyk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    sumYkTerms = Sykyk - Chat*Sxkyk - Sxkyk'*Chat' + Chat*Sxkxk*Chat';      
                    Sxkxkm1 = Sxkm1xk';
                    Sykxk = Sxkyk';

                    

                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    sumYkTerms=0.5*(sumYkTerms+sumYkTerms');
                    if(mPPCOEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(mPPCOEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    else
                        ScoreAMc=[];
                    end

                    ScorC=Rhat\(Sykxk-Chat*Sxkxk);
                    ScoreCMc=reshape(ScorC',numel(ScorC),1);

                    if(mPPCOEM_Constraints.QhatDiag)
                        if(mPPCOEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end


                    ScoreAlphaMc = sum(Rhat\(y-Chat*x_K-alphahat*ones(1,N)),2);
                    if(mPPCOEM_Constraints.RhatDiag)
                        if(mPPCOEM_Constraints.RhatIsotropic)
                            ScoreR  =-.5*(K*Dy*Rhat(1,1)^(-1) - Rhat(1,1)^(-2)*trace(sumYkTerms));
                        else
                            ScoreR  =(-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat)));
                        end
                        ScoreRMc = diag(ScoreR);
                    else
                        ScoreR   =-.5*(Rhat\(K*eye(size(Rhat)) - sumYkTerms/Rhat));
                        ScoreRMc =reshape(ScoreR',numel(ScoreR),1);
                    end


                    if(mPPCOEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);
                    ScoreMuMc=zeros(numCells,1);
                    ScoreBetaMc=[];
                    ScoreGammaMc=[];
                    % Cell Scores
                    for nc=1:numCells
                        if(strcmp(fitType,'poisson'))
                            Hk = (HkAll(:,:,nc));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms);
                            ScoreMuMc(nc) = sum(dN(nc,:)-ld,2);
                            ScoreBetaMc = [ScoreBetaMc; sum(repmat((dN(nc,:)-ld),[Dx 1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-ld,[nHist 1]).*Hk',2)];
                        elseif(strcmp(fitType,'binomial'))
                            Hk = (HkAll(:,:,nc));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
                            nHist = size(Hk,2);
                            
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ScoreMuMc(nc) = sum(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,2);
                            ScoreBetaMc = [ScoreBetaMc;sum(repmat(dN(nc,:).*(1-ld) - ld.*(1-ld),[Dx,1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,[nHist 1]).*Hk',2)];
                        end
                        
                    end
                    ScoreVec = [ScoreAMc; ScoreQMc; ScoreCMc; ScoreRMc];
                    if(mPPCOEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(mPPCOEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreAlphaMc];
                    ScoreVec = [ScoreVec; ScoreMuMc; ScoreBetaMc];
                    if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                        ScoreVec=[ScoreVec;ScoreGammaMc];
                    end
                    
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end
            end
            IMissing = 1/Mc*sum(IMc,3);
            IObs  = IComp-IMissing;  
            invIObs = eye(size(IObs))/IObs;
%             figure(1); subplot(1,2,1); imagesc(invIObs); subplot(1,2,2); imagesc(nearestSPD(invIObs));
            invIObs = nearestSPD(invIObs); % Find the nearest positive semidefinite approximation for the variance matrix
            VarVec = (diag(invIObs));
            SEVec = sqrt(VarVec);
            SEAterms = SEVec(1:n1);
            SEQterms = SEVec(n1+1:(n1+n2));
            SECterms = SEVec(n1+n2+1:(n1+n2+n3));
            SERterms = SEVec(n1+n2+n3+1:(n1+n2+n3+n4));
            SEPx0terms=SEVec(n1+n2+n3+n4+1:(n1+n2+n3+n4+n5));
            SEx0terms=SEVec(n1+n2+n3+n4+n5+1:(n1+n2+n3+n4+n5+n6));
            SEAlphaterms=SEVec(n1+n2+n3+n4+n5+n6+1:(n1+n2+n3+n4+n5+n6+n7));
            SEMuTerms = SEVec(n1+n2+n3+n4+n5+n6+n7+1:(n1+n2+n3+n4+n5+n6+n7+n8));
            SEBetaTerms = SEVec(n1+n2+n3+n4+n5+n6+n7+n8+1:(n1+n2+n3+n4+n5+n6+n7+n8+n9)); 
            SEGammaTerms = SEVec(n1+n2+n3+n4+n5+n6+n7+n8+n9+1:(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10)); 
            if(mPPCOEM_Constraints.EstimatePx0==1)
                SES = diag(SEPx0terms);
            end
            if(mPPCOEM_Constraints.Estimatex0==1)
                SEx0=SEx0terms;
            end

            if(mPPCOEM_Constraints.EstimateA==1)
                if(mPPCOEM_Constraints.AhatDiag==1)
                    SEA=diag(SEAterms);
                else
                    SEA=reshape(SEAterms,size(Ahat,2),size(Ahat,1))';
                end
            end
            SEC=reshape(SECterms,size(Chat,2),size(Chat,1))';
            SEAlpha=reshape(SEAlphaterms,size(alphahat,2),size(alphahat,1))';

            if(mPPCOEM_Constraints.RhatDiag==1)
                SER=diag(SERterms);
            else
                SER=reshape(SERterms,size(Rhat,2),size(Rhat,1))';
            end
            if(mPPCOEM_Constraints.QhatDiag==1)
                SEQ=diag(SEQterms);
            else
                SEQ=reshape(SEQterms,size(Qhat,2),size(Qhat,1))'; 
            end
            if(mPPCOEM_Constraints.EstimateA==1)
                SE.A = SEA;
            end
            SE.Q = SEQ;
            SE.C = SEC;
            SE.R = SER;
            SE.alpha = SEAlpha;

            if(mPPCOEM_Constraints.EstimatePx0==1)
                SE.Px0=SES;
            end
            if(mPPCOEM_Constraints.Estimatex0==1)
                SE.x0=SEx0;
            end
            
            SEMu = SEMuTerms;
            SEBeta=reshape(SEBetaTerms,size(betahat,2),size(betahat,1))';

            SE.mu = SEMu;
            SE.beta = SEBeta;
            if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                SEGamma=reshape(SEGammaTerms,size(gammahat,2),size(gammahat,1))';
                SE.gamma = SEGamma;
            end
            % Compute parameter p-values
            if(mPPCOEM_Constraints.EstimateA==1)
                clear h p;
                if(mPPCOEM_Constraints.AhatDiag==1)
                    VecParams = diag(Ahat);
                    VecSE     = diag(SEA);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pA = diag(p);
                else
                    VecParams = reshape(Ahat,[numel(Ahat) 1]);
                    VecSE     = reshape(SEA, [numel(Ahat) 1]);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end  
                    pA = reshape(p, [size(Ahat,1) size(Ahat,2)]);
                end
            end

            %C matrix
            clear h p;
            VecParams = reshape(Chat,[numel(Chat) 1]);
            VecSE     = reshape(SEC, [numel(Chat) 1]);
            for i=1:length(VecParams)
               [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end 
            pC = reshape(p, [size(Chat,1) size(Chat,2)]);

            %R matrix
            clear h p;
            if(mPPCOEM_Constraints.RhatDiag==1)
                if(mPPCOEM_Constraints.RhatIsotropic==1)
                    VecParams = Rhat(1,1);
                    VecSE     = SER(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pR = diag(p);
                else
                    VecParams = diag(Rhat);
                    VecSE     = diag(SER);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pR = diag(p);
                end
            else
                VecParams = reshape(Rhat,[numel(Rhat) 1]);
                VecSE     = reshape(SER, [numel(Rhat) 1]);
                for i=1:length(VecParams)
                   [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end  
                pR = reshape(p, [size(Rhat,1) size(Rhat,2)]);
            end

            %Q matrix
            clear h p;
            if(mPPCOEM_Constraints.QhatDiag==1)
                if(mPPCOEM_Constraints.QhatIsotropic==1)
                    VecParams = Qhat(1,1);
                    VecSE     = SEQ(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pQ = diag(p);
                else
                    VecParams = diag(Qhat);
                    VecSE     = diag(SEQ);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pQ = diag(p);
                end
            else
                VecParams = reshape(Qhat,[numel(Qhat) 1]);
                VecSE     = reshape(SEQ, [numel(Qhat) 1]);
                for i=1:length(VecParams)
                   [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end  
                pQ = reshape(p, [size(Qhat,1) size(Qhat,2)]);
            end
            %Px0
            if(mPPCOEM_Constraints.EstimatePx0==1)
                clear h p;
                if(mPPCOEM_Constraints.Px0Isotropic==1)
                    VecParams = Px0hat(1,1);
                    VecSE     = SES(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pPX0 = diag(p);
                else
                    VecParams = diag(Px0hat);
                    VecSE     = diag(SES);
                    for i=1:length(VecParams)
                        [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pPX0 = diag(p);
                end
            end

            clear h p;
            VecParams = alphahat;
            VecSE     = SEAlpha;
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end
            pAlpha = p';

            if(mPPCOEM_Constraints.Estimatex0==1)
                clear h p;
                VecParams = x0hat;
                VecSE     = SEx0;
                for i=1:length(VecParams)
                    [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end
                pX0 = p';
            end
            
            %Mu
            clear h p;
            VecParams = muhat;
            VecSE     = SEMu;
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end
            pMu = p';
            
            %Beta
            clear h p;
            VecParams = reshape(betahat,[numel(betahat),1]);
            VecSE     = reshape(SEBeta, [numel(SEBeta),1]);
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end    
            pBeta = reshape(p, [size(betahat,1) size(betahat,2)]);
            
            %Gamma
            clear h p;
            if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                VecParams = reshape(gammahat,[numel(gammahat),1]);
                VecSE     = reshape(SEGamma, [numel(gammahat),1]);
                for i=1:length(VecParams)
                    [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end    
                pGamma = reshape(p, [size(gammahat,1) size(gammahat,2)]);
            end
            if(mPPCOEM_Constraints.EstimateA==1)
                Pvals.A = pA;
            end
            Pvals.Q = pQ;
            Pvals.C = pC;
            Pvals.R = pR;
            Pvals.alpha = pAlpha;
            if(mPPCOEM_Constraints.EstimatePx0==1)
                Pvals.Px0 = pPX0;
            end
            if(mPPCOEM_Constraints.Estimatex0==1)
                Pvals.x0 = pX0;
            end
            Pvals.mu = pMu;
            Pvals.beta = pBeta;
            
            if(numel(gammahat)==1)
                if(gammahat~=0)
                    Pvals.gamma = pGamma;
                end
            else
                Pvals.gamma = pGamma;
            end

        end
        function [xKFinal,WKFinal,Ahat, Qhat, Chat, Rhat,alphahat, muhat, betahat, gammahat, x0hat, Px0hat, IC, SE, Pvals]=mPPCO_EM(y,dN, Ahat0, Qhat0, Chat0, Rhat0, alphahat0, mu, beta, fitType,delta, gamma, windowTimes, x0, Px0,mPPCOEM_Constraints,MstepMethod)
            numStates = size(Ahat0,1);
            if(nargin<17 || isempty(MstepMethod))
               MstepMethod='GLM'; %or NewtonRaphson 
            end
            if(nargin<16 || isempty(mPPCOEM_Constraints))
                mPPCOEM_Constraints = DecodingAlgorithms.mPPCO_EMCreateConstraints;
            end
            if(nargin<15 || isempty(Px0))
                Px0=10e-10*eye(numStates,numStates);
            end
            if(nargin<14 || isempty(x0))
                x0=zeros(numStates,1);
            end
            
            if(nargin<13 || isempty(windowTimes))
                if(isempty(gamma))
                    windowTimes =[];
                else
    %                 numWindows =length(gamma0)+1; 
                    windowTimes = 0:delta:(length(gamma)+1)*delta;
                end
            end
            if(nargin<12)
                gamma=[];
            end
            if(nargin<11 || isempty(delta))
                delta = .001;
            end
            if(nargin<10)
                fitType = 'poisson';
            end
            
            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
%                     HkAll{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                    HkAll(:,:,k) = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
%                     HkAll{k} = 0;
                    HkAll(:,:,k) = 0;
                end
                gamma=0;
            end



    %         tol = 1e-3; %absolute change;
            tolAbs = 1e-3;
            tolRel = 1e-3;
            llTol  = 1e-3;
            cnt=1;

            maxIter = 100;

            
            A0 = Ahat0;
            Q0 = Qhat0;
            C0 = Chat0;
            R0 = Rhat0;
            alpha0 = alphahat0;
           
            Ahat{1} = A0;
            Qhat{1} = Q0;
            Chat{1} = C0;
            Rhat{1} = R0;
            x0hat{1} = x0;
            Px0hat{1} = Px0;
            alphahat{1} = alpha0;
            muhat{1} = mu;
            betahat{1} = beta;
            gammahat{1} = gamma;
            yOrig=y;
            numToKeep=10;
            scaledSystem=1;
            
            if(scaledSystem==1)
                Tq = eye(size(Qhat{1}))/(chol(Qhat{1}));
                Tr = eye(size(Rhat{1}))/(chol(Rhat{1}));
                Ahat{1}= Tq*Ahat{1}/Tq;
                Chat{1}= Tr*Chat{1}/Tq;
                Qhat{1}= Tq*Qhat{1}*Tq';
                Rhat{1}= Tr*Rhat{1}*Tr';
                y= Tr*y;
                x0hat{1} = Tq*x0;
                Px0hat{1} = Tq*Px0*Tq';
                alphahat{1}= Tr*alphahat{1};  
                betahat{1}=(betahat{1}'/Tq)';
            end

            cnt=1;
            dLikelihood(1)=inf;
%             x0hat = x0;
            negLL=0;
            IkedaAcc=mPPCOEM_Constraints.EnableIkeda;
            %Forward EM
            stoppingCriteria =0;
%             logllNew= -inf;

            disp('                        Joint Point-Process/Gaussian Observation EM Algorithm                        ');     
            while(stoppingCriteria~=1 && cnt<=maxIter)
                 storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
                 storeIndP1= mod(cnt,numToKeep)+1;
                 storeIndM1= mod(cnt-2,numToKeep)+1;
                disp('--------------------------------------------------------------------------------------------------------');
                disp(['Iteration #' num2str(cnt)]);
                disp('--------------------------------------------------------------------------------------------------------');
                
                
                [x_K{storeInd},W_K{storeInd},ll(cnt),ExpectationSums{storeInd}]=...
                    DecodingAlgorithms.mPPCO_EStep(Ahat{storeInd},Qhat{storeInd},Chat{storeInd},Rhat{storeInd}, y, alphahat{storeInd},dN, muhat{storeInd}, betahat{storeInd},fitType,delta,gammahat{storeInd},HkAll, x0hat{storeInd}, Px0hat{storeInd});
                
                [Ahat{storeIndP1}, Qhat{storeIndP1}, Chat{storeIndP1}, Rhat{storeIndP1}, alphahat{storeIndP1}, muhat{storeIndP1}, betahat{storeIndP1}, gammahat{storeIndP1},x0hat{storeIndP1},Px0hat{storeIndP1}] ...
                    = DecodingAlgorithms.mPPCO_MStep(dN, y,x_K{storeInd},W_K{storeInd},x0hat{storeInd}, Px0hat{storeInd}, ExpectationSums{storeInd}, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,mPPCOEM_Constraints,MstepMethod);
              
                if(IkedaAcc==1)
                    disp(['****Ikeda Acceleration Step****']);
                    %y=Cx+alpha+wk wk~Normal with covariance Rk
                     ykNew = mvnrnd((Chat{storeIndP1}*x_K{storeInd}+alphahat{storeIndP1}*ones(1,size(x_K{storeInd},2)))',Rhat{storeIndP1})';
                     
%                      if(gammahat{storeIndP1}==0)% No history effect
%                         dataMat = [ones(size(y,2),1) x_K{storeInd}']; % design matrix: X 
%                         coeffsMat = [muhat{storeIndP1} betahat{storeIndP1}']; % coefficient vector: beta
%                         minTime=0;
%                         maxTime=(size(dN,2)-1)*delta;
%                         time=minTime:delta:maxTime;
%                         clear nstNew;
%                         for cc=1:length(muhat{storeIndP1})
%                              tempData  = exp(dataMat*coeffsMat(cc,:)');
% 
%                              if(strcmp(fitType,'poisson'))
%                                  lambdaData = tempData;
%                              else
%                                 lambdaData = tempData./(1+tempData); % Conditional Intensity Function for ith cell
%                              end
%                              lambda{cc}=Covariate(time,lambdaData./delta, ...
%                                  '\Lambda(t)','time','s','spikes/sec',...
%                                  {strcat('\lambda_{',num2str(cc),'}')},{{' ''b'' '}});
%                              lambda{cc}=lambda{cc}.resample(1/delta);
% 
%                              % generate one realization for each cell
%                              tempSpikeColl{cc} = CIF.simulateCIFByThinningFromLambda(lambda{cc},1);          
%                              nstNew{cc} = tempSpikeColl{cc}.getNST(1);     % grab the realization
%                              nstNew{cc}.setName(num2str(cc));              % give each cell a unique name
% %                              subplot(4,3,[8 11]);
% %                              h2=lambda{cc}.plot([],{{' ''k'', ''LineWidth'' ,.5'}}); 
% %                              legend off; hold all; % Plot the CIF
% 
%                         end
%                         
%                         spikeColl = nstColl(nstNew); % Create a neural spike train collection
%                      else
%                          time;
%                      end
                     
                     dNNew=dN;%spikeColl.dataToMatrix';
                     %dNNew(dNNew>1)=1; % more than one spike per bin will be treated as one spike. In
                                    % general we should pick delta small enough so that there is
                                    % only one spike per bin
                                    
                                    
                                    
                     [x_KNew,W_KNew,llNew,ExpectationSumsNew]=...
                        DecodingAlgorithms.mPPCO_EStep(Ahat{storeInd},Qhat{storeInd},Chat{storeInd},Rhat{storeInd}, ykNew, alphahat{storeInd},dNNew, muhat{storeInd}, betahat{storeInd},fitType,delta,gammahat{storeInd},HkAll, x0, Px0);

                
                     [AhatNew, QhatNew, ChatNew, RhatNew, alphahatNew, muhatNew, betahatNew, gammahatNew,x0new,Px0new] ...
                        = DecodingAlgorithms.mPPCO_MStep(dNNew, ykNew,x_KNew,W_KNew, x0hat{storeInd}, Px0hat{storeInd}, ExpectationSumsNew, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,mPPCOEM_Constraints,MstepMethod);
               
                    Ahat{storeIndP1} = 2*Ahat{storeIndP1}-AhatNew;
                    Qhat{storeIndP1} = 2*Qhat{storeIndP1}-QhatNew;
                    Qhat{storeIndP1} = (Qhat{storeIndP1}+Qhat{storeIndP1}')/2;
                    Chat{storeIndP1} = 2*Chat{storeIndP1}-ChatNew;
                    Rhat{storeIndP1} = 2*Rhat{storeIndP1}-RhatNew;
                    Rhat{storeIndP1} = (Rhat{storeIndP1}+Rhat{storeIndP1}')/2;
                    alphahat{storeIndP1}=2*alphahat{storeIndP1}-alphahatNew;
%                     muhat{storeIndP1}= 2*muhat{storeIndP1}-muhatNew;
%                     betahat{storeIndP1} = 2*betahat{storeIndP1}-betahatNew;
%                     gammahat{storeIndP1}= 2*gammahat{storeIndP1}-gammahatNew;
%                     x0hat{storeIndP1}   = 2*x0hat{storeIndP1} - x0new;
%                     Px0hat{storeIndP1}  = 2*Px0hat{storeIndP1}- Px0new;
%                     [V,D] = eig(Px0hat{storeIndP1});
%                     D(D<0)=1e-9;
%                     Px0hat{storeIndP1} = V*D*V';
%                     Px0hat{storeIndP1}  = (Px0hat{storeIndP1}+Px0hat{storeIndP1}')/2;
                    
               
                end
                if(mPPCOEM_Constraints.EstimateA==0)
                    Ahat{storeIndP1}=Ahat{storeInd};
                end
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(ll(cnt)-ll(cnt-1));%./abs(ll(cnt-1));
                end
                if(cnt==1)
                    QhatInit = Qhat{1};
                    RhatInit = Rhat{1};
                    xKInit = x_K{1};
                end
                %Plot the progress
%                 if(mod(cnt,2)==0)
                if(cnt==1)
                    scrsz = get(0,'ScreenSize');
                    h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
                end
                    figure(h);
                    time = linspace(minTime,maxTime,size(x_K{storeInd},2));
                    subplot(2,5,[1 2 6 7]); plot(1:cnt,ll,'k','Linewidth', 2); hy=ylabel('Log Likelihood'); hx=xlabel('Iteration'); axis auto;
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    subplot(2,5,3:5); hNew=plot(time, x_K{storeInd}','Linewidth', 2); hy=ylabel('States'); hx=xlabel('time [s]');
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(time, xKInit','--','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    subplot(2,5,8); hNew=plot(diag(Qhat{storeInd}),'o','Linewidth', 2); hy=ylabel('Q'); hx=xlabel('Diagonal Entry');
                    set(gca, 'XTick'       , 1:1:length(diag(Qhat{storeInd})));
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(diag(QhatInit),'r.','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    subplot(2,5,9); hNew=plot(diag(Rhat{storeInd}),'o','Linewidth', 2); hy=ylabel('R'); hx=xlabel('Diagonal Entry');
                    set(gca, 'XTick'       , 1:1:length(diag(Rhat{storeInd})));
                    set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    hold on; hOrig=plot(diag(RhatInit),'r.','Linewidth', 2); 
                    legend([hOrig(1) hNew(1)],'Initial','Current');
                    
                    
                    subplot(2,5,10); imagesc(Rhat{storeInd}); ht=title('R Matrix Image'); 
                    set(gca, 'XTick'       , 1:1:length(diag(Rhat{storeInd})), 'YTick', 1:1:length(diag(Rhat{storeInd})));
                    set(ht,'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                    drawnow;
                    hold off;
%                 end
                
                if(cnt==1)
                    dMax=inf;
                else
                 dQvals = max(max(abs(sqrt(Qhat{storeInd})-sqrt(Qhat{storeIndM1}))));
                 dRvals = max(max(abs(sqrt(Rhat{storeInd})-sqrt(Rhat{storeIndM1}))));
                 dAvals = max(max(abs((Ahat{storeInd})-(Ahat{storeIndM1}))));
                 dCvals = max(max(abs((Chat{storeInd})-(Chat{storeIndM1}))));
                 dMuvals = max(abs((muhat{storeInd})-(muhat{storeIndM1})));
                 dAlphavals = max(abs((alphahat{storeInd})-(alphahat{storeIndM1})));
                 dBetavals = max(max(abs((betahat{storeInd})-(betahat{storeIndM1}))));
                 dGammavals = max(max(abs((gammahat{storeInd})-(gammahat{storeIndM1}))));
                 dMax = max([dQvals,dRvals,dAvals,dCvals,dMuvals,dAlphavals,dBetavals,dGammavals]);
                end

% 
%                 dQRel = max(abs(dQvals./sqrt(Qhat(:,storeIndM1))));
%                 dGammaRel = max(abs(dGamma./gammahat(storeIndM1,:)));
%                 dMaxRel = max([dQRel,dGammaRel]);
                if(cnt==1)
                    disp(['Max Parameter Change: N/A']);
                else
                    disp(['Max Parameter Change: ' num2str(dMax)]);
                end
                cnt=(cnt+1);
                if(dMax<tolAbs)
                    stoppingCriteria=1;
                    display(['         EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
                    negLL=0;
                end
            
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
                    stoppingCriteria=1;
                    display(['         EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
                    
                    negLL=1;
                end
                

            end
            
            disp('--------------------------------------------------------------------------------------------------------');


            maxLLIndex  = find(ll == max(ll),1,'first');
            maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
            if(maxLLIndex==1)
%                 maxLLIndex=cnt-1;
                maxLLIndex =1;
                maxLLIndMod = 1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
               maxLLIndMod = 1;
%             else
%                maxLLIndMod = mod(maxLLIndex,numToKeep); 
               
            end
            nIter   = cnt-1;  
%             maxLLIndMod
           
            xKFinal = x_K{maxLLIndMod};
            WKFinal = W_K{maxLLIndMod};
            Ahat = Ahat{maxLLIndMod};
            Qhat = Qhat{maxLLIndMod};
            Chat = Chat{maxLLIndMod};
            Rhat = Rhat{maxLLIndMod};
            alphahat = alphahat{maxLLIndMod};
            muhat= muhat{maxLLIndMod};
            betahat = betahat{maxLLIndMod};
            gammahat = gammahat{maxLLIndMod};
            x0hat =x0hat{maxLLIndMod};
            Px0hat=Px0hat{maxLLIndMod};
            
             if(scaledSystem==1)
               Tq = eye(size(Qhat))/(chol(Q0));
               Tr = eye(size(Rhat))/(chol(R0));
               Ahat=Tq\Ahat*Tq;
               Qhat=(Tq\Qhat)/Tq';
               Chat=Tr\Chat*Tq;
               Rhat=(Tr\Rhat)/Tr';
               alphahat=Tr\alphahat;
               xKFinal = Tq\xKFinal;
               x0hat = Tq\x0hat;
               Px0hat= (Tq\Px0hat)/(Tq');
               tempWK =zeros(size(WKFinal));
               for kk=1:size(WKFinal,3)
                tempWK(:,:,kk)=(Tq\WKFinal(:,:,kk))/Tq';
               end
               WKFinal = tempWK;
               betahat=(betahat'*Tq)';
             end
            llFinal=ll(end);
            ll = ll(maxLLIndex);
            ExpectationSumsFinal = ExpectationSums{maxLLIndMod};
% AhatNew, QhatNew, ChatNew, RhatNew, alphahatNew, muhatNew, betahatNew, gammahatNew,x0new,Px0new
            if(nargout>13)
                [SE, Pvals]=DecodingAlgorithms.mPPCO_ComputeParamStandardErrors(y, dN,...
                    xKFinal, WKFinal, Ahat, Qhat, Chat, Rhat, alphahat, x0hat, Px0hat, ExpectationSumsFinal,...
                    fitType, muhat, betahat, gammahat, windowTimes, HkAll,...
                    mPPCOEM_Constraints);
            end
            
            %Compute number of parameters
            if(mPPCOEM_Constraints.EstimateA==1 && mPPCOEM_Constraints.AhatDiag==1)
                n1=size(Ahat,1); 
            elseif(mPPCOEM_Constraints.EstimateA==1 && mPPCOEM_Constraints.AhatDiag==0)
                n1=numel(Ahat);
            else 
                n1=0;
            end
            if(mPPCOEM_Constraints.QhatDiag==1 && mPPCOEM_Constraints.QhatIsotropic==1)
                n2=1;
            elseif(mPPCOEM_Constraints.QhatDiag==1 && mPPCOEM_Constraints.QhatIsotropic==0)
                n2=size(Qhat,1);
            else
                n2=numel(Qhat);
            end

            n3=numel(Chat); 
            if(mPPCOEM_Constraints.RhatDiag==1 && mPPCOEM_Constraints.RhatIsotropic==1)
                n4=1;
            elseif(mPPCOEM_Constraints.QhatDiag==1 && mPPCOEM_Constraints.QhatIsotropic==0)
                n4=size(Rhat,1);
            else
                n4=numel(Rhat);
            end

            if(mPPCOEM_Constraints.EstimatePx0==1 && mPPCOEM_Constraints.Px0Isotropic==1)
                n5=1;
            elseif(mPPCOEM_Constraints.EstimatePx0==1 && mPPCOEM_Constraints.Px0Isotropic==0)
                n5=size(Px0hat,1);
            else
                n5=0;
            end

            if(mPPCOEM_Constraints.Estimatex0==1)   
                n6=size(x0hat,1);
            else
                n6=0;
            end

            n7=size(alphahat,1);
            n8=size(muhat,1);
            n9=numel(betahat);
            if(numel(gammahat)==1)
                if(gammahat==0)
                    n10=0;
                else
                    n10=1;
                end
            else
                n10=numel(gammahat);
            end
            nTerms=n1+n2+n3+n4+n5+n6+n7+n8+n9+n10;
            
            K  = size(y,2); 
            Dx = size(Ahat,2);
            sumXkTerms = ExpectationSums{maxLLIndMod}.sumXkTerms;
            llobs = ll + Dx*K/2*log(2*pi)+K/2*log(det(Qhat))...
                + 1/2*trace(Qhat\sumXkTerms)...
                + Dx/2*log(2*pi)+1/2*log(det(Px0hat)) ...
                + 1/2*Dx;
            AIC = 2*nTerms - 2*llobs;
            AICc= AIC+ 2*nTerms*(nTerms+1)/(K-nTerms-1);
            BIC = -2*llobs+nTerms*log(K);
            IC.AIC = AIC;
            IC.AICc= AICc;
            IC.BIC = BIC;
            IC.llobs = llobs;
            IC.llcomp=ll;
         
            
        end
        function [x_K,W_K,logll,ExpectationSums]=mPPCO_EStep(A,Q,C,R, y, alpha,dN, mu, beta,fitType,delta,gamma,HkAll, x0, Px0)
             DEBUG = 0;

             minTime=0;
             maxTime=(size(dN,2)-1)*delta;


   
            [numCells,K]   = size(dN); 
            Dx = size(A,2);
            Dy = size(C,1);
            x_p     = zeros( size(A,2), K );
            x_u     = zeros( size(A,2), K );
            W_p    = zeros( size(A,2),size(A,2), K);
            W_u    = zeros( size(A,2),size(A,2), K );
            

            [x_p, W_p, x_u, W_u] = DecodingAlgorithms.mPPCODecodeLinear(A, Q, C, R, y, alpha, dN,mu,beta,fitType,delta,gamma,[],x0,Px0,HkAll);
            
            [x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u);
            
            %Best estimates of initial states given the data
            W1G0 = A*Px0*A' + Q;
            L0=Px0*A'/W1G0;
            
            Ex0Gy = x0+L0*(x_K(:,1)-x_p(:,1));        
            Px0Gy = Px0+L0*(eye(size(W_K(:,:,1)))/(W_K(:,:,1))-eye(size(W1G0))/W1G0)*L0';
            Px0Gy = (Px0Gy+Px0Gy')/2;
            numStates = size(x_K,1);
            Wku=zeros(numStates,numStates,K,K);
            Tk = zeros(numStates,numStates,K-1);
            for k=1:K
                Wku(:,:,k,k)=W_K(:,:,k);
            end

            for u=K:-1:2
                for k=(u-1):-1:(u-1)
                    Tk(:,:,k)=A;
%                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
                    Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
                    Wku(:,:,u,k)=Wku(:,:,k,u)';
                end
            end
            
            %All terms
            Sxkm1xk = zeros(Dx,Dx);
            Sxkxkm1 = zeros(Dx,Dx);
            Sxkm1xkm1 = zeros(Dx,Dx);
            Sxkxk = zeros(Dx,Dx);
            Sykyk = zeros(Dy,Dy);
            Sxkyk = zeros(Dx,Dy);
            for k=1:K
                if(k==1)
                    Sxkm1xk   = Sxkm1xk+Px0*A'/W_p(:,:,1)*Wku(:,:,1,1);
                    Sxkm1xkm1 = Sxkm1xkm1+Px0+x0*x0';     
                else
%                   
                      Sxkm1xk =  Sxkm1xk+Wku(:,:,k-1,k)+x_K(:,k-1)*x_K(:,k)';
                       
                      Sxkm1xkm1= Sxkm1xkm1+Wku(:,:,k-1,k-1)+x_K(:,k-1)*x_K(:,k-1)';
                end
                Sxkxk = Sxkxk+Wku(:,:,k,k)+x_K(:,k)*x_K(:,k)';
                Sykyk = Sykyk+(y(:,k)-alpha)*(y(:,k)-alpha)';
                Sxkyk = Sxkyk+x_K(:,k)*(y(:,k)-alpha)';

            end
            Sx0x0 = Px0+x0*x0';
            Sxkxk = 0.5*(Sxkxk+Sxkxk');
            Sykyk = 0.5*(Sykyk+Sykyk');
            sumXkTerms = Sxkxk-A*Sxkm1xk-Sxkm1xk'*A'+A*Sxkm1xkm1*A';
            sumYkTerms = Sykyk - C*Sxkyk - Sxkyk'*C' + C*Sxkxk*C';      
            Sxkxkm1 = Sxkm1xk';
            
%             if(strcmp(fitType,'poisson'))
%                 sumPPll=0;
%                 for c=1:numCells
%                     Hk=HkAll{c};
%                     for k=1:K
%                         xk = x_K(:,k);
%                         if(numel(gamma)==1)
%                             gammaC=gamma;
%                         else 
%                             gammaC=gamma(:,c);
%                         end
%                         terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
%                         Wk = W_K(:,:,k);
%                         ld = exp(terms);
%                         bt = beta(:,c);
%                         ExplambdaDelta =ld+0.5*trace(bt*bt'*ld*Wk);
%                         ExplogLD = terms;
%                         sumPPll=sumPPll+dN(c,k).*ExplogLD - ExplambdaDelta;
%                     end
%                   
%                             
%                 end
%             elseif(strcmp(fitType,'binomial'))
%                 sumPPll=0;
%                 for c=1:numCells
%                     Hk=HkAll{c};
%                     for k=1:K
%                         xk = x_K(:,k);
%                         if(numel(gamma)==1)
%                             gammaC=gamma;
%                         else 
%                             gammaC=gamma(:,c);
%                         end
%                         terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
%                         Wk = W_K(:,:,k);
%                         ld = exp(terms)./(1+exp(terms));
%                         bt = beta(:,c);
%                         ExplambdaDelta =ld+0.5*trace(bt*bt'*ld*(1-ld)*(1-2*ld)*Wk);
%                         ExplogLD = log(ld)+0.5*trace(-(bt*bt'*ld*(1-ld))*Wk);
%                         sumPPll=sumPPll+dN(c,k).*ExplogLD - ExplambdaDelta;
%                     end
%                   
%                             
%                 end
%             end
            %Vectorize for loop over cells
            if(strcmp(fitType,'poisson'))
                sumPPll=0;
                HkPerm =permute(HkAll,[2 3 1]);
                for k=1:K
%                    Hk=squeeze(HkAll(k,:,:)); 
                   Hk = HkPerm(:,:,k);
                   if(size(Hk,1)==numCells)
                       Hk = Hk';
                   end
                   xk = x_K(:,k);
                   if(numel(gamma)==1)
                        gammaC=repmat(gamma,1,numCells);
                   else 
                        gammaC=gamma;
                   end
                   if(size(gammaC,2)~=numCells)
                       gammaC = repmat(gammaC,[1 numCells]);
                   end
                   terms=mu+beta'*xk+diag(gammaC'*Hk);
                   Wk = W_K(:,:,k);
                   ld = exp(terms);
                   bt = beta;
                   ExplambdaDelta =ld+0.5*(ld.*diag((bt'*Wk*bt)));
                   ExplogLD = terms;
                   sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta);
                        
                end
                
            %Vectorize over number of cells
            elseif(strcmp(fitType,'binomial'))
                sumPPll=0;
                HkPerm = permute(HkAll,[2 3 1]);
                for k=1:K
%                     Hk=squeeze(HkAll(k,:,:)); 
                    HkPerm = HkPerm(:,:,k);
                    if(size(Hk,1)==numCells)
                       Hk = Hk';
                    end
                    xk = x_K(:,k);
                    if(numel(gamma)==1)
                        gammaC=repmat(gamma,1,numCells);
                    else 
                        gammaC=gamma;
                    end
                    if(size(gammaC,2)~=numCells)
                       gammaC = repmat(gammaC,[1 numCells]);
                    end
                   terms=mu+beta'*xk+diag(gammaC'*Hk);
                   Wk = W_K(:,:,k);
                   ld = exp(terms)./(1+exp(terms));
                   bt = beta;     
                   ExplambdaDelta = ld+0.5*(ld.*(1-ld).*(1-2.*ld)).*diag((bt'*Wk*bt));
                   ExplogLD = log(ld)+0.5*(-ld.*(1-ld)).*diag(bt'*Wk*bt);
                   sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta); 
                    
                end

                
            end

            logll = -Dx*K/2*log(2*pi)-K/2*log(det(Q))-Dy*K/2*log(2*pi) ...
                    -K/2*log(det(R))- Dx/2*log(2*pi) -1/2*log(det(Px0))  ...
                    +sumPPll - 1/2*trace((eye(size(Q))/Q)*sumXkTerms) ...
                    -1/2*trace((eye(size(R))/R)*sumYkTerms) ...
                    -Dx/2;
                string0 = ['logll: ' num2str(logll)];
                disp(string0);
                if(DEBUG==1)
                    string1 = ['-K/2*log(det(Q)):' num2str(-K/2*log(det(Q)))];
                    string11 = ['-K/2*log(det(R)):' num2str(-K/2*log(det(R)))];
                    string12= ['Constants: ' num2str(-Dx*K/2*log(2*pi)-Dy*K/2*log(2*pi)- Dx/2*log(2*pi) -Dx/2 -1/2*log(det(Px0)))];
                    string2 = ['SumPPll: ' num2str(sumPPll)];
                    string3 = ['-.5*trace(Q\sumXkTerms): ' num2str(-.5*trace(Q\sumXkTerms))];
                    string4 = ['-.5*trace(R\sumYkTerms): ' num2str(-.5*trace(R\sumYkTerms))];

                    disp(string1);
                    disp(['Q=' num2str(diag(Q)')]);
                    disp(string11);
                    disp(['R=' num2str(diag(R)')]);
                    disp(string12);
                    disp(string2);
                    disp(string3);
                    disp(string4);
                end

                ExpectationSums.Sxkm1xkm1=Sxkm1xkm1;
                ExpectationSums.Sxkm1xk=Sxkm1xk;
                ExpectationSums.Sxkxkm1=Sxkxkm1;
                ExpectationSums.Sxkxk=Sxkxk;
                ExpectationSums.Sxkyk=Sxkyk;
                ExpectationSums.Sykyk=Sykyk;
                ExpectationSums.sumXkTerms=sumXkTerms;
                ExpectationSums.sumYkTerms=sumYkTerms;
                ExpectationSums.sumPPll=sumPPll;
                ExpectationSums.Sx0 = Ex0Gy;
                ExpectationSums.Sx0x0 = Px0Gy + Ex0Gy*Ex0Gy';

        end
        function [Ahat, Qhat, Chat, Rhat, alphahat, muhat_new, betahat_new, gammahat_new, x0hat, Px0hat] = mPPCO_MStep(dN, y,x_K,W_K,x0, Px0, ExpectationSums,fitType, muhat, betahat,gammahat, windowTimes, HkAll,mPPCOEM_Constraints,MstepMethod)
            if(nargin<14 || isempty(MstepMethod))
                MstepMethod = 'GLM'; %GLM or NewtonRaphson
            end
            if(nargin<13 || isempty(mPPCOEM_Constraints))
                mPPCOEM_Constraints = DecodingAlgorithms.mPPCO_EMCreateConstraints;
            end
           
            Sxkm1xkm1=ExpectationSums.Sxkm1xkm1;
            Sxkm1xk=ExpectationSums.Sxkm1xk;
            Sxkxkm1=ExpectationSums.Sxkxkm1;
            Sxkxk=ExpectationSums.Sxkxk;
            Sxkyk=ExpectationSums.Sxkyk;
            Sykyk=ExpectationSums.Sykyk;
            sumXkTerms = ExpectationSums.sumXkTerms;
            sumYkTerms = ExpectationSums.sumYkTerms;
            Sx0 = ExpectationSums.Sx0;
            Sx0x0 = ExpectationSums.Sx0x0;
            [dx,K] = size(x_K);   
            dy=size(y,1);
            numCells=size(dN,1);
            
            if(mPPCOEM_Constraints.AhatDiag==1)
                I=eye(dx,dx);
                Ahat = (Sxkxkm1.*I)/(Sxkm1xkm1.*I);
            else
                Ahat = Sxkxkm1/Sxkm1xkm1;
            end
            Chat = Sxkyk'/Sxkxk;             
            alphahat = sum(y - Chat*x_K,2)/K;
            
            if(mPPCOEM_Constraints.QhatDiag==1)
                 if(mPPCOEM_Constraints.QhatIsotropic==1)
                     Qhat=1/(dx*K)*trace(sumXkTerms)*eye(dx,dx);
                 else
                     I=eye(dx,dx);
                     Qhat=1/K*(sumXkTerms.*I);
                     Qhat = (Qhat + Qhat')/2;
                 end
             else
                 Qhat=1/K*sumXkTerms;
                 Qhat = (Qhat + Qhat')/2;
             end
             dy=size(sumYkTerms,1);
             if(mPPCOEM_Constraints.RhatDiag==1)
                 if(mPPCOEM_Constraints.RhatIsotropic==1)
                     I=eye(dy,dy);
                     Rhat = 1/(dy*K)*trace(sumYkTerms)*I;
                 else
                     
                     I=eye(dy,dy);
                     Rhat = 1/K*(sumYkTerms.*I);
                     Rhat = (Rhat + Rhat')/2;
                 end
             else
                 Rhat = 1/K*(sumYkTerms);
                 Rhat = (Rhat + Rhat')/2;  
             end
             if(mPPCOEM_Constraints.Estimatex0)
                x0hat = (inv(Px0)+Ahat'/Qhat*Ahat)\(Ahat'/Qhat*x_K(:,1)+Px0\x0);
            else
                x0hat = x0;
            end
             
            if(mPPCOEM_Constraints.EstimatePx0==1)
                if(mPPCOEM_Constraints.Px0Isotropic==1)
                   Px0hat=(trace(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0'))/(dx*K))*eye(dx,dx); 
                else
                    I=eye(dx,dx);
                    Px0hat =(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0')).*I;
                    Px0hat = (Px0hat+Px0hat')/2;
                end
                
            else
                Px0hat =Px0;
            end
             
             betahat_new =betahat;
             gammahat_new = gammahat;
             muhat_new = muhat;
             
            %Compute the new CIF beta using the GLM
            if(strcmp(fitType,'poisson'))
                algorithm = 'GLM';
            else
                algorithm = 'BNLRCG';
            end
            
            % Estimate params via GLM

            if(strcmp(MstepMethod,'GLM'))
                clear c; close all;
                time=(0:length(x_K)-1)*.001;
                labels = cell(1,dx);
                labels2 = cell(1,dx+1);
                labels2{1} = 'vel';
                for i=1:dx
                    labels{i} = strcat('v',num2str(i));
                    labels2{i+1} = strcat('v',num2str(i));
                end
                vel = Covariate(time,x_K','vel','time','s','m/s',labels);
                baseline = Covariate(time,ones(length(time),1),'Baseline','time','s','',...
                    {'constant'});
                for i=1:size(dN,1)
                    spikeTimes = time(find(dN(i,:)==1));
                    nst{i} = nspikeTrain(spikeTimes);
                end
                nspikeColl = nstColl(nst);
                cc = CovColl({vel,baseline});
                trial = Trial(nspikeColl,cc);
                selfHist = windowTimes ; NeighborHist = []; sampleRate = 1000; 
                clear c;
                
                

                if(gammahat==0)
                    c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,[],NeighborHist); 
                else
                    c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,selfHist,NeighborHist); 
                end
                c{1}.setName('Baseline');
                cfgColl= ConfigColl(c);
                warning('OFF');

                results = Analysis.RunAnalysisForAllNeurons(trial,cfgColl,0,algorithm);
                temp = FitResSummary(results);
                tempCoeffs = squeeze(temp.getCoeffs);
                if(gammahat==0)
                    betahat(1:dx,:) = tempCoeffs(2:(dx+1),:);
                    muhat = tempCoeffs(1,:)';
                else
                    betahat(1:dx,:) = tempCoeffs(2:(dx+1),:);
                    muhat = tempCoeffs(1,:)';
                    histTemp = squeeze(temp.getHistCoeffs);
                    histTemp = reshape(histTemp, [length(windowTimes)-1 numCells]);
                    histTemp(isnan(histTemp))=0;
                    gammahat=histTemp;
                end
            else
                
                
            % Estimate via Newton-Raphson
                 fprintf(['****M-step for beta**** \n']);
                 McExp=50;    
                 xKDrawExp = zeros(size(x_K,1),K,McExp);
                 diffTol = 1e-5;

                % Generate the Monte Carlo samples
                for k=1:K
                    WuTemp=(W_K(:,:,k));
                    [chol_m,p]=chol(WuTemp);
                    z=normrnd(0,1,size(x_K,1),McExp);
                    xKDrawExp(:,k,:)=repmat(x_K(:,k),[1 McExp])+(chol_m*z);
                end
                % Stimulus Coefficients
                pool = matlabpool('size');
                if(pool==0)
                    xkPerm = permute(xKDrawExp,[1 3 2]);
                    for c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
                        fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
                            if(iter==1)
                                fprintf('%d',iter);
                            else
                                fprintf(',%d',iter);
                            end
                            if(strcmp(fitType,'poisson'))
                                HessianTerm = zeros(size(x_K,1),size(x_K,1));
                                GradTerm = zeros(size(x_K,1),1);
                                xkPerm = permute(xKDrawExp,[1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaXkXkT = 1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    GradTerm = GradTerm+dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                    HessianTerm=HessianTerm-ExpLambdaXkXkT;

                                end

                            elseif(strcmp(fitType,'binomial'))
                                HessianTerm = zeros(size(x_K,1),size(x_K,1));
                                GradTerm = zeros(size(x_K,1),1);
                                xkPerm = permute(xKDrawExp,[1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaSquaredXk = 1/McExp*sum(repmat(ld.^2,[size(xk,1),1]).*xk,2);
                                    GradTerm = GradTerm+dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                    HessianTerm=HessianTerm+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;

                                end

                            end
                            if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                betahat_newTemp = betahat_new(:,c);
                            else
                                betahat_newTemp = (betahat_new(:,c)-HessianTerm\GradTerm);
                                if(any(isnan(betahat_newTemp)))
                                    betahat_newTemp = betahat_new(:,c);

                                end
                            end
                            mabsDiff = max(abs(betahat_newTemp - betahat_new(:,c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            betahat_new(:,c)=betahat_newTemp;
                            iter=iter+1;
                        end
                        fprintf('\n');              
                    end 
                else
                    HessianTerm = zeros(size(betahat,1),size(betahat,1),numCells);
                    GradTerm = zeros(size(betahat,1),numCells);
                    betahat_newTemp=betahat_new;
                    xkPerm = permute(xKDrawExp,[1 3 2]);
                    parfor c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
                        fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
                            if(iter==1)
                                fprintf('%d',iter);
                            else
                                fprintf(',%d',iter);
                            end

                            if(strcmp(fitType,'poisson'))

                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaXkXkT = 1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    if(k==1)
                                        GradTerm(:,c) = dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                        HessianTerm(:,:,c)=-ExpLambdaXkXkT;
                                    else
                                        GradTerm(:,c) = GradTerm(:,c)+dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                        HessianTerm(:,:,c)=HessianTerm(:,:,c)-ExpLambdaXkXkT;
                                    end

                                end

                            elseif(strcmp(fitType,'binomial'))

                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaSquaredXk = 1/McExp*sum(repmat(ld.^2,[size(xk,1),1]).*xk,2);
                                    if(k==1)
                                        GradTerm(:,c) = dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                        HessianTerm(:,:,c)=ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;
                                    else
                                        GradTerm(:,c) = GradTerm(:,c)+dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                        HessianTerm(:,:,c)=HessianTerm(:,:,c)+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;
                                    end
                                end

                            end
                            if(any(any(isnan(HessianTerm(:,:,c)))) || any(any(isinf(HessianTerm(:,:,c)))))
                                betahat_newTemp = betahat_new(:,c);
                            else
                                betahat_newTemp = (betahat_new(:,c)-HessianTerm(:,:,c)\GradTerm(:,c));
                                if(any(isnan(betahat_newTemp)))
                                    betahat_newTemp = betahat_new(:,c);

                                end
                            end
                            mabsDiff = max(abs(betahat_newTemp - betahat_new(:,c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            betahat_new(:,c)=betahat_newTemp;
                            iter=iter+1;
                        end
                        fprintf('\n');              
                    end 
                end
                clear GradTerm HessianTerm;
                 %Compute the CIF means 
                 if(pool==0)
                     xkPerm = permute(xKDrawExp,[1 3 2]);
                     for c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
    %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
    %                         if(iter==1)
    %                             fprintf('%d',iter);
    %                         else
    %                             fprintf(',%d',iter);
    %                         end
                            if(strcmp(fitType,'poisson'))
                                HessianTerm = zeros(size(1,1),size(1,1));
                                GradTerm = zeros(size(1,1),1);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaDelta = 1/McExp*sum(ld,2);
                                    GradTerm = GradTerm+(dN(c,k) - ExpLambdaDelta);
                                    HessianTerm=HessianTerm-ExpLambdaDelta;

                                end

                            elseif(strcmp(fitType,'binomial'))
                                HessianTerm = zeros(size(1,1),size(1,1));
                                GradTerm = zeros(size(1,1),1);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExpLambdaDelta =1/McExp*(sum(ld,2));
                                    ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                    ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                    GradTerm = GradTerm+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                    HessianTerm=HessianTerm+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);

                                end

                            end
                            if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                muhat_newTemp = muhat_new(c);
                            else
                                muhat_newTemp = (muhat_new(c)-HessianTerm\GradTerm);
                                if(any(isnan(muhat_newTemp)))
                                    muhat_newTemp = muhat_new(c);

                                end
                            end
                            mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            muhat_new(c)=muhat_newTemp;
                            iter=iter+1;
                        end
    %                     fprintf('\n');              
                     end 
                 else
                    HessianTerm = zeros(1,numCells);
                    GradTerm = zeros(1,numCells);
                    xkPerm = permute(xKDrawExp,[1 3 2]);
                    parfor c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
    %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
    %                         if(iter==1)
    %                             fprintf('%d',iter);
    %                         else
    %                             fprintf(',%d',iter);
    %                         end
                            if(strcmp(fitType,'poisson'))
                                for k=1:K
                                    Hk = squeeze(HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaDelta = 1/McExp*sum(ld,2);
                                    if(k==1)
                                        GradTerm(c) = (dN(c,k) - ExpLambdaDelta);
                                        HessianTerm(c)=-ExpLambdaDelta;
                                    else
                                        GradTerm(c) = GradTerm(c)+(dN(c,k) - ExpLambdaDelta);
                                        HessianTerm(c)=HessianTerm(c)-ExpLambdaDelta;
                                    end

                                end

                            elseif(strcmp(fitType,'binomial'))
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExpLambdaDelta =1/McExp*(sum(ld,2));
                                    ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                    ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                    if(k==1)
                                        GradTerm(c) = (dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                        HessianTerm(c)=(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);
                                    else
                                         GradTerm(c) = GradTerm(c)+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                        HessianTerm(c)=HessianTerm(c)+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);
                                    end

                                end

                            end
                            if(any(any(isnan(HessianTerm(c)))) || any(any(isinf(HessianTerm(c)))))
                                muhat_newTemp = muhat_new(c);
                            else
                                muhat_newTemp = (muhat_new(c)-HessianTerm(c)\GradTerm(c));
                                if(any(isnan(muhat_newTemp)))
                                    muhat_newTemp = muhat_new(c);

                                end
                            end
                            mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            muhat_new(c)=muhat_newTemp;
                            iter=iter+1;
                        end
    %                     fprintf('\n');              
                     end 
                 end
                 clear HessianTerm GradTerm;
                 
                 
                 %Compute the history coeffs
                 if(~isempty(windowTimes) && any(any(gammahat_new~=0)))
                     if(pool==0)
                         xkPerm = permute(xKDrawExp,[1 3 2]);
                         for c=1:numCells
                            converged=0;
                            iter = 1;
                            maxIter=100;
        %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                            while(~converged && iter<maxIter)
        %                         if(iter==1)
        %                             fprintf('%d',iter);
        %                         else
        %                             fprintf(',%d',iter);
        %                         end
                            
                                if(strcmp(fitType,'poisson'))
                                    HessianTerm = zeros(size(gammahat,1),size(gammahat,1));
                                    GradTerm = zeros(size(gammahat,1),1);
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end

                                        if(numel(gammahat_new)==1)
                                            gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                                        else 
                                            gammaC=gammahat_new(:,c);
                                        end

                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms);
                                        ExpLambdaDelta = 1/McExp*sum(ld,2);
                                        GradTerm = GradTerm+(dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                        HessianTerm=HessianTerm-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);

                                    end

                                elseif(strcmp(fitType,'binomial'))
                                    HessianTerm = zeros(size(gammahat,1),size(gammahat,1));
                                    GradTerm = zeros(size(gammahat,1),1);
                                    for k=1:K
                                        Hk = squeeze(HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end

                                        if(numel(gammahat_new)==1)
                                            gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                                        else 
                                            gammaC=gammahat_new(:,c);
                                        end

                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms)./(1+exp(terms));
                                        ExpLambdaDelta =1/McExp*(sum(ld,2));
                                        ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                        ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                        GradTerm = GradTerm+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                        HessianTerm=HessianTerm+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);

                                    end

                                end
                                if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                    gammahat_newTemp = gammahat_new(:,c);
                                else
                                    gammahat_newTemp = (gammahat_new(:,c)-HessianTerm\GradTerm);
                                    if(any(isnan(gammahat_newTemp)))
                                        gammahat_newTemp = gammahat_new(:,c);

                                    end
                                end
                                mabsDiff = max(abs(gammahat_newTemp - gammahat_new(:,c)));
                                if(mabsDiff<diffTol)
                                    converged=1;
                                end
                                gammahat_new(:,c)=gammahat_newTemp;
                                iter=iter+1;
                            end
        %                     fprintf('\n');              
                         end 
                     else
                         HessianTerm = zeros(size(gammahat,1),size(gammahat,1),numCells);
                         GradTerm = zeros(size(gammahat,1),numCells);
                         xkPerm = permute(xKDrawExp,[1 3 2]);
                         parfor c=1:numCells
                            converged=0;
                            iter = 1;
                            maxIter=100;
        %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                            if(numel(gammahat_new)==1)
                                gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat_new(:,c);
                            end
                            
                            while(~converged && iter<maxIter)
        %                         if(iter==1)
        %                             fprintf('%d',iter);
        %                         else
        %                             fprintf(',%d',iter);
        %                         end

                                if(strcmp(fitType,'poisson'))
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end



                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms);
                                        ExpLambdaDelta = 1/McExp*sum(ld,2);
                                        if(k==1)
                                            GradTerm(:,c) = (dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                            HessianTerm(:,:,c)=-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);
                                        else
                                            GradTerm(:,c) = GradTerm(:,c)+(dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                            HessianTerm(:,:,c)=HessianTerm(:,:,c)-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);
                                        end
                                    end

                                elseif(strcmp(fitType,'binomial'))
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end


                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms)./(1+exp(terms));
                                        ExpLambdaDelta =1/McExp*(sum(ld,2));
                                        ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                        ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                        if(k==1)
                                            GradTerm(:,c) = (dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                            HessianTerm(:,:,c)=(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);
                                        else
                                            GradTerm(:,c) = GradTerm(:,c)+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                            HessianTerm(:,:,c)=HessianTerm(:,:,c)+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);
                                        end

                                    end

                                end
                                if(any(any(isnan(HessianTerm(:,:,c)))) || any(any(isinf(HessianTerm(:,:,c)))))
                                    gammahat_newTemp = gammaC;
                                else
                                    gammahat_newTemp = (gammaC-HessianTerm(:,:,c)\GradTerm(:,c));
                                    if(any(isnan(gammahat_newTemp)))
                                        gammahat_newTemp = gammaC;

                                    end
                                end
                                mabsDiff = max(abs(gammahat_newTemp - gammaC));
                                if(mabsDiff<diffTol)
                                    converged=1;
                                end
                                gammaC=gammahat_newTemp;
                                iter=iter+1;
                            end
                            gamma_new(:,c) =gammaC;
        %                     fprintf('\n');              
                         end 
                     end
                 end
                 clear HessianTerm GradTerm;       
                 
%                  muhat_new =muhat;
%                  for c=1:numCells
%                      converged=0;
%                      iter = 1;
%                      maxIter=100;
%                      while(~converged && iter<maxIter)
%                         if(strcmp(fitType,'poisson'))
%                             gradQ=zeros(size(muhat_new(c),2),1);
%                             jacQ =zeros(size(muhat_new(c),2),size(muhat_new(c),2));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                 ld = exp(terms);
%                                 bt = betahat(:,c);
%                                 ExplambdaDelta =ld +0.5*trace(ld*bt*bt'*Wk);
% 
% 
%                                 gradQ = gradQ + dN(c,k)' - ExplambdaDelta;
%                                 jacQ  = jacQ  - ExplambdaDelta;
%                             end
% 
% 
%                         elseif(strcmp(fitType,'binomial'))
%                             gradQ=zeros(size(muhat_new(c),2),1);
%                             jacQ =zeros(size(muhat_new(c),2),size(muhat_new(c),2));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                 ld = exp(terms)./(1+exp(terms));
%                                 bt = betahat(:,c);
%                                 ExplambdaDelta = ld+0.5*trace(bt*bt'*(ld)*(1-ld)*(1-2*ld)*Wk);
%                                 ExplambdaDeltaSq = (ld)^2+...
%                                     0.5*trace((ld)^2*(1-ld)*(2-3*ld)*bt*bt'*Wk);
%                                 ExplambdaDeltaCubed = (ld)^3+...
%                                     0.5*trace(3*(ld)^3*(3-7*ld+4*(ld)^2)*bt*bt'*Wk);
% 
%                                 gradQ = gradQ + dN(c,k)' -(dN(c,k)+1)*ExplambdaDelta...
%                                     +ExplambdaDeltaSq;
%                                 jacQ  = jacQ  - (dN(c,k)+1)*ExplambdaDelta...
%                                     +(dN(c,k)+3)*ExplambdaDeltaSq...
%                                     -3*ExplambdaDeltaCubed;
%                             end
% 
%                         end
%     %                     gradQ=0.01*gradQ;
%                         muhat_newTemp = (muhat_new(c)'-(1/jacQ)*gradQ)';
%                         if(any(isnan(muhat_newTemp)))
%                             muhat_newTemp = muhat_new(c);
% 
%                         end
%                         mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
%                         if(mabsDiff<10^-2)
%                             converged=1;
%                         end
%                         muhat_new(c)=muhat_newTemp;
%                         iter=iter+1;
%                      end
% 
%                 end
% 
%     %             Compute the history parameters
%                 gammahat_new = gammahat;
%                 if(~isempty(windowTimes) && any(any(gammahat_new~=0)))
%                      for c=1:numCells
%                          converged=0;
%                          iter = 1;
%                          maxIter=100;
%                          while(~converged && iter<maxIter)
%                             if(strcmp(fitType,'poisson'))
%                                 gradQ=zeros(size(gammahat_new(c),2),1);
%                                 jacQ =zeros(size(gammahat_new(c),2),size(gammahat_new(c),2));
%                                 for k=1:K
% %                                     Hk=HkAll{c};
%                                     Hk = squeeze(HkAll(:,:,c));
%                                     Wk = W_K(:,:,k);
%                                     if(numel(gammahat)==1)
%                                         gammaC=gammahat;
%                                     else 
%                                         gammaC=gammahat(:,c);
%                                     end
%                                     terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                     ld = exp(terms);
%                                     bt = betahat(:,c);
%                                     ExplambdaDelta =ld +0.5*trace(bt*bt'*ld*Wk);
% 
% 
%                                     gradQ = gradQ + (dN(c,k)' - ExplambdaDelta)*Hk;
%                                     jacQ  = jacQ  - ExplambdaDelta*Hk*Hk';
%                                 end
% 
% 
%                             elseif(strcmp(fitType,'binomial'))
%                                 gradQ=zeros(size(gammahat_new(c),2),1);
%                                 jacQ =zeros(size(gammahat_new(c),2),size(gammahat_new(c),2));
%                                 for k=1:K
% %                                     Hk=HkAll{c};
%                                     Hk = squeeze(HkAll(:,:,c));
%                                     Wk = W_K(:,:,k);
%                                     if(numel(gammahat)==1)
%                                         gammaC=gammahat;
%                                     else 
%                                         gammaC=gammahat(:,c);
%                                     end
%                                     terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                     ld = exp(terms)./(1+exp(terms));
%                                     bt = betahat(:,c);
%                                     ExplambdaDelta =ld...
%                                         +0.5*trace(bt*bt'*ld*(1-ld)*(1-2*ld)*Wk);
%                                     ExplambdaDeltaSq=ld^2 ...
%                                         +trace((ld^2*(1-ld)*(2-3*ld)*bt*bt')*Wk);
%                                     ExplambdaDeltaCubed=ld^3 ...
%                                         +0.5*trace((9*(ld^3)*(1-ld)^2*bt*bt'-3*(ld^4)*(1-ld)*bt*bt')*Wk);
%                                     gradQ = gradQ + (dN(c,k) - (dN(c,k)+1)*ExplambdaDelta+ExplambdaDeltaSq)*Hk;
%                                     jacQ  = jacQ  + -ExplambdaDelta*(dN(c,k)+1)*Hk*Hk'...
%                                         +ExplambdaDeltaSq*(dN(c,k)+3)*Hk*Hk'...
%                                         -ExplambdaDeltaCubed*2*Hk*Hk';
%                                 end
% 
%                             end
% 
% 
%     %                         gradQ=0.01*gradQ;
% 
%                             gammahat_newTemp = (gammahat_new(:,c)-(eye(size(Hk,2),size(Hk,2))/jacQ)*gradQ');
%                             if(any(isnan(gammahat_newTemp)))
%                                 gammahat_newTemp = gammahat_new(:,c);
% 
%                             end
%                             mabsDiff = max(abs(gammahat_newTemp - gammahat_new(:,c)));
%                             if(mabsDiff<10^-2)
%                                 converged=1;
%                             end
%                             gammahat_new(:,c)=gammahat_newTemp;
%                             iter=iter+1;
%                          end
% 
%                     end
%     %                  gammahat(:,c) = gammahat_new;
%                 end
%              betahat =betahat_new;
%              gammahat = gammahat_new;
%              muhat = muhat_new;
            end           
        end
    
        %% Point Process EM
        function C = PP_EMCreateConstraints(EstimateA, AhatDiag,QhatDiag,QhatIsotropic,Estimatex0,EstimatePx0, Px0Isotropic,mcIter, EnableIkeda)
            %By default, all parameters are estimated. To empose diagonal
            %structure on the EM parameter results must pass in the
            %constraints element
            if(nargin<9 || isempty(EnableIkeda))
                EnableIkeda=0;
            end
            if(nargin<8 || isempty(mcIter))
                mcIter=1000;
            end
            if(nargin<7 || isempty(Px0Isotropic))
                Px0Isotropic=0;
            end
            if(nargin<6 || isempty(EstimatePx0))
                EstimatePx0=1;
            end
            if(nargin<5 || isempty(Estimatex0))
                Estimatex0=1;
            end
            if(nargin<4 || isempty(QhatIsotropic))
                QhatIsotropic=0;
            end
            if(nargin<3 || isempty(QhatDiag))
                QhatDiag=1;
            end
            if(nargin<2)
                AhatDiag=0;
            end
            if(nargin<1)
                EstimateA=1;
            end
            C.EstimateA = EstimateA;
            C.AhatDiag = AhatDiag;
            C.QhatDiag = QhatDiag;
            if(QhatDiag && QhatIsotropic)
                C.QhatIsotropic=1;
            else
                C.QhatIsotropic=0;
            end
            C.Estimatex0 = Estimatex0;
            C.EstimatePx0 = EstimatePx0;
            if(EstimatePx0 && Px0Isotropic)
                C.Px0Isotropic=1;
            else
                C.Px0Isotropic=0; 
            end
            C.mcIter = mcIter;
            C.EnableIkeda=EnableIkeda;
        end  
        function [SE,Pvals,nTerms] = PP_ComputeParamStandardErrors(dN, xKFinal, WKFinal, Ahat, Qhat, x0hat, Px0hat, ExpectationSumsFinal, fitType, muhat, betahat, gammahat, windowTimes, HkAll, PPEM_Constraints)

            % Use inverse observed information matrix to estimate the standard errors of the estimated model parameters
         % Requires computation of the complete information matrix and an estimate of the missing information matrix

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Complete Information Matrices     
        % Recall from McLachlan and Krishnan Eq. 4.7
        %    Io(theta;y) = Ic(theta;y) - Im(theta;y)
        %    Io(theta;y) = Ic(theta;y) - cov(Sc(X;theta)Sc(X;theta)')
        % where Sc(X;theta) is the score vector of the complete log likelihood
        % function evaluated at theta. We first compute Ic term by term and then
        % approximate the covariance term using Monte Carlo approximation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(nargin<19 || isempty(PPEM_Constraints))
                PPEM_Constraints=DecodingAlgorithms.PP_EMCreateConstraints;
            end

            
            if(PPEM_Constraints.EstimateA==1)
                if(PPEM_Constraints.AhatDiag==1)
                    IAComp=zeros(numel(diag(Ahat)),numel(diag(Ahat)));
                else
                    IAComp=zeros(numel(Ahat),numel(Ahat));
                end
                [n1,n2] =size(Ahat);
                el=(eye(n1,n1));
                em=(eye(n2,n2));
                cnt=1;
                N=size(xKFinal,2);

                if(PPEM_Constraints.AhatDiag==1)
                    for l=1:n1
                        for m=l
                            termMat=Qhat\el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1.*eye(n1,n2);
                            termvec = diag(termMat);
                            IAComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                else
                    for l=1:n1
                        for m=1:n2
                            termMat=(inv(Qhat))*el(:,l)*em(:,m)'*ExpectationSumsFinal.Sxkm1xkm1;
                            termvec=reshape(termMat',1,numel(Ahat));
                            IAComp(:,cnt)=termvec';
                            cnt=cnt+1;
                        end
                    end
                end
            end

           
            [n1,n2] =size(Qhat);
            el=(eye(n1,n1));
            em=(eye(n2,n2));
            cnt=1;
            if(PPEM_Constraints.QhatDiag==1)
                if(PPEM_Constraints.QhatIsotropic==1)
                    IQComp=zeros(1,1);
                    IQComp =  0.5*N*dx*Qhat(1,1)^(-2); 
                else
                    IQComp=zeros(numel(diag(Qhat)),numel(diag(Qhat)));
                    for l=1:n1
                        for m=l
                            termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                            termvec=diag(termMat);
                            IQComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            else
                IQComp=zeros(numel(Qhat),numel(Qhat));
                for l=1:n1
                    for m=1:n2
                        termMat= N/2*(Qhat)\em(:,m)*el(:,l)'/(Qhat);
                        termvec=reshape(termMat',1,numel(Qhat));
                        IQComp(:,cnt)=termvec;
                        cnt=cnt+1;
                    end
                end
            end

            if(PPEM_Constraints.EstimatePx0==1)
                if(PPEM_Constraints.Px0Isotropic==1)
                    ISComp =  0.5*dx*Px0hat(1,1)^(-2);
                else
                    ISComp=zeros(numel(diag(Px0hat)),numel(diag(Px0hat)));
                    [n1,n2] =size(Px0hat);
                    el=(eye(n1,n1));
                    em=(eye(n2,n2));
                    cnt=1;
                    for l=1:n1
                        for m=l
                            termMat= 1/2*(Px0hat)\em(:,m)*el(:,l)'/(Px0hat);
                            termvec=diag(termMat);
                            ISComp(:,cnt)=termvec;
                            cnt=cnt+1;
                        end
                    end
                end
            end

            if(PPEM_Constraints.Estimatex0==1)
                Ix0Comp=eye(size(Px0hat))/Px0hat+(Ahat'/Qhat)*Ahat;
            end

            
            K=size(xKFinal,2);
            numCells=size(betahat,2);
            McExp=PPEM_Constraints.mcIter; 
            xKDrawExp = zeros(size(xKFinal,1),K,McExp);
            

            % Generate the Monte Carlo
            for k=1:K
                WuTemp=squeeze(WKFinal(:,:,k));
                [chol_m,p]=chol(WuTemp);
                z=normrnd(0,1,size(xKFinal,1),McExp);
                xKDrawExp(:,k,:)=repmat(xKFinal(:,k),[1 McExp])+(chol_m*z);
            end
            
            IBetaComp =zeros(size(xKFinal,1)*numCells,size(xKFinal,1)*numCells);
            xkPerm = permute(xKDrawExp,[1 3 2]);
            pools = matlabpool('size'); %number of parallel workers 
           if(strcmp(fitType,'poisson'))
                for c=1:numCells
                    HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1),K);
                    parfor k=1:K
%                         Hk = squeeze(HkAll(:,:,c));
                        Hk = (HkAll(k,:,c));
                        Wk = WKFinal(:,:,k);
                       
%                         xk = squeeze(xKDrawExp(:,k,:));
                        xk=xkPerm(:,:,k);
                       if(size(Hk,1)==numCells)
                           Hk = Hk';
                       end
                   
                        if(numel(gammahat)==1)
                            gammaC=gammahat;
%                             gammaC=repmat(gammaC,[1 numCells]);
                        else 
                            gammaC=gammahat(:,c);
                        end

                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                        ld=exp(terms);
                        
                        HessianTerm(:,:,k)=-1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                    end
                    startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                    IBetaComp(startInd:endInd,startInd:endInd)=-sum(HessianTerm,3);
                end
            else
                for c=1:numCells
                    HessianTerm = zeros(size(xKFinal,1),size(xKFinal,1),K);
                    parfor k=1:K
%                         Hk = squeeze(HkAll(:,:,c));
                        Hk = (HkAll(k,:,c));
                        Wk = WKFinal(:,:,k);
%                         xk = squeeze(xKDrawExp(:,k,:));
                        xk = (xkPerm(:,:,k));
                        if(size(Hk,1)==numCells)
                           Hk = Hk';
                        end
                   
                        if(numel(gammahat)==1)
                            gammaC=gammahat;
%                             gammaC=repmat(gammaC,[1 numCells]);
                        else 
                            gammaC=gammahat(:,c);
                        end
                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                        ld=exp(terms)./(1+exp(terms));
                        ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                        ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                        ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                        HessianTerm(:,:,k)+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;
                        
                    end
                    startInd = size(betahat,1)*(c-1)+1; endInd = size(betahat,1)*c;
                    IBetaComp(startInd:endInd,startInd:endInd)=-sum(HessianTerm,3);
                end
            end

            
            %CIF means
            IMuComp=zeros(numel(muhat),numel(muhat));
            xkPerm = permute(xKDrawExp,[1 3 2]);
            if(pools==0)
                for c=1:numCells
                    if(strcmp(fitType,'poisson'))
                        HessianTerm = 0;
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(:,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                            ld = exp(terms);
                            HessianTerm=HessianTerm-1/McExp*sum(ld,2);
                        end
                    elseif(strcmp(fitType,'binomial'))
                        HessianTerm = 0;
                        for k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(:,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                            ld = exp(terms)./(1+exp(terms));
                            ExplambdaDelta = 1/McExp*sum(ld,2);
                            ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                            ExplambdaDeltaCubed = 1/McExp*sum(ld.^3,2);
                            HessianTerm = HessianTerm -(dN(c,k)+1)*ExplambdaDelta ...
                                +(dN(c,k)+3)*ExplambdaDeltaSquare-3*ExplambdaDeltaCubed;
                        end
                    end
                    IMuComp(c,c) = -HessianTerm;
                end
            else
                for c=1:numCells
                    if(strcmp(fitType,'poisson'))
                        HessianTerm = zeros(K,1);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld = exp(terms);
                            HessianTerm(k)=-1/McExp*sum(ld,2);
                        end
                    elseif(strcmp(fitType,'binomial'))
                        HessianTerm = zeros(K,1);
                        parfor k=1:K
    %                         Hk = squeeze(HkAll(:,:,c));
                            Hk = (HkAll(k,:,c));
                            if(size(Hk,1)==numCells)
                               Hk = Hk';
                            end
    %                         xk = squeeze(xKDrawExp(:,k,:));
                            xk = xkPerm(:,:,k);
                            Wk = WKFinal(:,:,k);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,c);
                            end
                            terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ExplambdaDelta = 1/McExp*sum(ld,2);
                            ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                            ExplambdaDeltaCubed = 1/McExp*sum(ld.^3,2);
                            HessianTerm(k) =  -(dN(c,k)+1)*ExplambdaDelta ...
                                +(dN(c,k)+3)*ExplambdaDeltaSquare-3*ExplambdaDeltaCubed;
                        end
                    end
                    IMuComp(c,c) = -sum(HessianTerm);
                end
            end
            
            
                       % Gamma Information Matrix
            IGammaComp = zeros(numel(gammahat),numel(gammahat));
            if(~isempty(windowTimes) && any(any(gammahat~=0)))
                xkPerm = permute(xKDrawExp,[1 3 2]);
                if(pools==0)
                     for c=1:numCells
                       if(strcmp(fitType,'poisson'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2));
                            for k=1:K
    %                             Hk = squeeze(HkAll(:,:,c));
                                Hk = (HkAll(:,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                ld = exp(terms);
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                HessianTerm=HessianTerm-Hk(k,:)'*Hk(k,:)*ExplambdaDelta;
                            end
                       elseif(strcmp(fitType,'binomial'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2));
                            for k=1:K
                                Hk = (HkAll(:,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                ld = exp(terms)./(1+exp(terms));
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                                ExplambdaDeltaCubed  = 1/McExp*sum(ld.^2,2);
                                HessianTerm=HessianTerm+(-ExplambdaDelta*(dN(c,k)+1)...
                                    +ExplambdaDeltaSquare*(dN(c,k)+3)...
                                    -2*ExplambdaDeltaCubed)*Hk(k,:)'*Hk(:,k);
                            end
                       end
                       startInd=size(HkAll,2)*(c-1)+1; endInd = size(HkAll,2)*c;
                       IGammaComp(startInd:endInd,startInd:endInd) = -HessianTerm;
                     end

                else
            
                    for c=1:numCells
                       if(strcmp(fitType,'poisson'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2),K);
                            parfor k=1:K
    %                             Hk = squeeze(HkAll(:,:,c));
                                Hk = (HkAll(k,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                                ld = exp(terms);
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                HessianTerm(:,:,k)=-Hk'*Hk*ExplambdaDelta;
                            end
                        elseif(strcmp(fitType,'binomial'))
                            HessianTerm = zeros(size(HkAll,2),size(HkAll,2),K);

                            parfor k=1:K
                                Hk = (HkAll(k,:,c));
                                if(size(Hk,1)==numCells)
                                   Hk = Hk';
                                end
    %                             xk = squeeze(xKDrawExp(:,k,:));
                                xk = xkPerm(:,:,k);
                                Wk = WKFinal(:,:,k);
                                if(numel(gammahat)==1)
                                    gammaC=gammahat;
                                else 
                                    gammaC=gammahat(:,c);
                                end
                                terms=muhat(c)+betahat(:,c)'*xk+gammaC'*Hk';
                                ld = exp(terms)./(1+exp(terms));
                                ExplambdaDelta = 1/McExp*sum(ld,2);
                                ExplambdaDeltaSquare = 1/McExp*sum(ld.^2,2);
                                ExplambdaDeltaCubed  = 1/McExp*sum(ld.^2,2);
                                HessianTerm(:,:,k)=+(-ExplambdaDelta*(dN(c,k)+1)...
                                    +ExplambdaDeltaSquare*(dN(c,k)+3)...
                                    -2*ExplambdaDeltaCubed)*Hk'*Hk;
                            end
                       end
                       startInd=size(HkAll,2)*(c-1)+1; endInd = size(HkAll,2)*c;
                       IGammaComp(startInd:endInd,startInd:endInd) = -sum(HessianTerm,3);
                    end

                end
            end
        
              
            
            if(PPEM_Constraints.EstimateA==1)
                n1=size(IAComp,1); 
            else
                n1=0;
            end
            n2=size(IQComp,1); 
          
            if(PPEM_Constraints.EstimatePx0==1)
                n3=size(ISComp,1); 
            else
                n3=0;
            end
            if(PPEM_Constraints.Estimatex0==1)   
                n4=size(Ix0Comp,1);
            else
                n4=0;
            end
            n5=size(IMuComp,1);
            n6=size(IBetaComp,1);
            if(numel(gammahat)==1)
                if(gammahat==0)
                    n7=0;
                end
            else
                n7=size(IGammaComp,1);
            end
            nTerms=n1+n2+n3+n4+n5+n6+n7;
            IComp = zeros(nTerms,nTerms);
            if(PPEM_Constraints.EstimateA==1)
                IComp(1:n1,1:n1)=IAComp;
            end
            offset=n1+1;
            IComp(offset:(n1+n2),offset:(n1+n2))=IQComp;
            offset=n1+n2+1;
            if(PPEM_Constraints.EstimatePx0==1);
                IComp(offset:(n1+n2+n3),offset:(n1+n2+n3))=ISComp;
            end
            offset=n1+n2+n3+1;
            if(PPEM_Constraints.Estimatex0==1)
                IComp(offset:(n1+n2+n3+n4),offset:(n1+n2+n3+n4))=Ix0Comp;
            end
            offset=n1+n2+n3+n4+1;
            IComp(offset:(n1+n2+n3+n4+n5),offset:(n1+n2+n3+n4+n5))=IMuComp;
            offset=n1+n2+n3+n4+n5+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6),offset:(n1+n2+n3+n4+n5+n6))=IBetaComp;
            offset=n1+n2+n3+n4+n5+n6+1;
            IComp(offset:(n1+n2+n3+n4+n5+n6+n7),offset:(n1+n2+n3+n4+n5+n6+n7))=IGammaComp;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Missing Information Matrix
            %Approximate cov(Sc(X;theta)Sc(X;theta)')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Mc=PPEM_Constraints.mcIter;
            xKDraw = zeros(size(xKFinal,1),N,Mc);

            % Generate the Monte Carlo samples for the unobserved data
            for n=1:N
                WuTemp=(WKFinal(:,:,n));
                [chol_m,p]=chol(WuTemp);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                xKDraw(:,n,:)=repmat(xKFinal(:,n),[1 Mc])+(chol_m*z);
            end


            if(PPEM_Constraints.EstimatePx0|| PPEM_Constraints.Estimatex0)
                [chol_m,p]=chol(Px0hat);
                z=normrnd(0,1,size(xKFinal,1),Mc);
                x0Draw=repmat(x0hat,[1 Mc])+(chol_m*z); 
            else
               x0Draw=repmat(x0hat, [1 Mc]);

            end

            IMc = zeros(nTerms,nTerms,Mc);
            % Emperically estimate the covariance of the score
            pools = matlabpool('size'); %number of parallel workers 
            if(pools==0) % parallel toolbox is not enabled;
                for c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                       
                    end
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    Sxkxkm1 = Sxkm1xk';
                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    if(PPEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(PPEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    else
                        ScoreAMc=[];
                    end

              
                    if(PPEM_Constraints.QhatDiag)
                        if(PPEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end

                    if(PPEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);
                    ScoreMuMc=zeros(numCells,1);
                    ScoreBetaMc=[];
                    ScoreGammaMc=[];
                    % Cell Scores
                    for nc=1:numCells
                        if(strcmp(fitType,'poisson'))
                            Hk = (HkAll(:,:,nc));
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms);
                            ScoreMuMc(nc) = sum(dN(nc,:)-ld,2);
                            ScoreBetaMc = [ScoreBetaMc; sum(repmat((dN(nc,:)-ld),[Dx 1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-ld,[nHist 1]).*Hk',2)];
                        elseif(strcmp(fitType,'binomial'))
                            Hk = (HkAll(:,:,nc));
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ScoreMuMc(nc) = sum(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,2);
                            ScoreBetaMc = [ScoreBetaMc;sum(repmat(dN(nc,:).*(1-ld) - ld.*(1-ld),[Dx,1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,[nHist 1]).*Hk',2)];
                        end
                        
                    end
                    ScoreVec = [ScoreAMc; ScoreQMc];
                    if(PPEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(PPEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreMuMc; ScoreBetaMc];
                    if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                            ScoreVec=[ScoreVec;ScoreGammaMc];
                    end
                    
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end
            else %Use the parallel toolbox
                parfor c=1:Mc
                    x_K=xKDraw(:,:,c);
                    x_0=x0Draw(:,c);

                    Dx=size(x_K,1);
                    Sxkm1xk = zeros(Dx,Dx);
                    Sxkm1xkm1 = zeros(Dx,Dx);
                    Sxkxk = zeros(Dx,Dx);

                    for k=1:K
                        if(k==1)
                            Sxkm1xk   = Sxkm1xk+x_0*x_K(:,k)';
                            Sxkm1xkm1 = Sxkm1xkm1+x_0*x_0';     
                        else
                            Sxkm1xk =  Sxkm1xk+x_K(:,k-1)*x_K(:,k)';
                            Sxkm1xkm1= Sxkm1xkm1+x_K(:,k-1)*x_K(:,k-1)';
                        end
                        Sxkxk = Sxkxk+x_K(:,k)*x_K(:,k)';
                       
                    end
                    Sxkxk = 0.5*(Sxkxk+Sxkxk');
                    sumXkTerms = Sxkxk-Ahat*Sxkm1xk-Sxkm1xk'*Ahat'+Ahat*Sxkm1xkm1*Ahat';
                    Sxkxkm1 = Sxkm1xk';
                    sumXkTerms=0.5*(sumXkTerms+sumXkTerms');
                    ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                    if(PPEM_Constraints.EstimateA==1)
                        ScorA=Qhat\(Sxkxkm1-Ahat*Sxkm1xkm1);
                        if(PPEM_Constraints.AhatDiag==1)
                            ScoreAMc=diag(ScorA);
                        else
                            ScoreAMc=reshape(ScorA',numel(Ahat),1);
                        end
                    else
                        ScoreAMc=[];
                    end


              
                    if(PPEM_Constraints.QhatDiag)
                        if(PPEM_Constraints.QhatIsotropic)
                            ScoreQ  =-.5*(K*Dx*Qhat(1,1)^(-1) - Qhat(1,1)^(-2)*trace(sumXkTerms));
                        else
                            ScoreQ  =(-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat)));
                        end
                        ScoreQMc = diag(ScoreQ);
                    else
                        ScoreQ   =-.5*(Qhat\(K*eye(size(Qhat)) - sumXkTerms/Qhat));
                        ScoreQMc =reshape(ScoreQ',numel(ScoreQ),1);
                    end

                    if(PPEM_Constraints.Px0Isotropic==1)
                        ScoreSMc=-.5*(Dx*Px0hat(1,1)^(-1) - Px0hat(1,1)^(-2)*trace((x_0-x0hat)*(x_0-x0hat)'));
                    else
                        ScorS  =-.5*(Px0hat\(eye(size(Px0hat)) - (x_0-x0hat)*(x_0-x0hat)'/Px0hat));
                        ScoreSMc = diag(ScorS);
                    end

                    Scorx0=(-Px0hat\(x_0-x0hat))+Ahat'/Qhat*(x_K(:,1)-Ahat*x_0);
                    Scorex0Mc=reshape(Scorx0',numel(Scorx0),1);
                    ScoreMuMc=zeros(numCells,1);
                    ScoreBetaMc=[];
                    ScoreGammaMc=[];
                    % Cell Scores
                    for nc=1:numCells
                        if(strcmp(fitType,'poisson'))
                            Hk = (HkAll(:,:,nc));
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms);
                            ScoreMuMc(nc) = sum(dN(nc,:)-ld,2);
                            ScoreBetaMc = [ScoreBetaMc; sum(repmat((dN(nc,:)-ld),[Dx 1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-ld,[nHist 1]).*Hk',2)];
                        elseif(strcmp(fitType,'binomial'))
                            Hk = (HkAll(:,:,nc));
                            nHist = size(Hk,2);
                            if(numel(gammahat)==1)
                                gammaC=gammahat;
                            else 
                                gammaC=gammahat(:,nc);
                            end
                            terms=muhat(nc)+betahat(:,nc)'*x_K+gammaC'*Hk';
                            ld = exp(terms)./(1+exp(terms));
                            ScoreMuMc(nc) = sum(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,2);
                            ScoreBetaMc = [ScoreBetaMc;sum(repmat(dN(nc,:).*(1-ld) - ld.*(1-ld),[Dx,1]).*x_K,2)];
                            ScoreGammaMc= [ScoreGammaMc;sum(repmat(dN(nc,:)-(dN(nc,:)+1).*ld+ld.^2,[nHist 1]).*Hk',2)];
                        end
                        
                    end
                    ScoreVec = [ScoreAMc; ScoreQMc];
                    if(PPEM_Constraints.EstimatePx0==1)
                        ScoreVec = [ScoreVec; ScoreSMc]; 
                    end
                    if(PPEM_Constraints.Estimatex0==1)
                        ScoreVec = [ScoreVec; Scorex0Mc];
                    end
                    ScoreVec = [ScoreVec; ScoreMuMc; ScoreBetaMc];
                    if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                            ScoreVec=[ScoreVec;ScoreGammaMc];
                    end
                    
                    IMc(:,:,c)=ScoreVec*ScoreVec';    
                end

            end
            IMissing = 1/Mc*sum(IMc,3);
            IObs  = IComp-IMissing;  
            invIObs = eye(size(IObs))/IObs;
%             figure(1); subplot(1,2,1); imagesc(invIObs); subplot(1,2,2); imagesc(nearestSPD(invIObs));
            invIObs = nearestSPD(invIObs); % Find the nearest positive semidefinite approximation for the variance matrix
            VarVec = (diag(invIObs));
            SEVec = sqrt(VarVec);
            SEAterms = SEVec(1:n1);
            SEQterms = SEVec(n1+1:(n1+n2));
            SEPx0terms=SEVec(n1+n2+1:(n1+n2+n3));
            SEx0terms=SEVec(n1+n2+n3+1:(n1+n2+n3+n4));
            SEMuTerms = SEVec(n1+n2+n3+n4+1:(n1+n2+n3+n4+n5));
            SEBetaTerms = SEVec(n1+n2+n3+n4+n5+1:(n1+n2+n3+n4+n5+n6)); 
            SEGammaTerms = SEVec(n1+n2+n3+n4+n5+n6+1:(n1+n2+n3+n4+n5+n6+n7)); 
            if(PPEM_Constraints.EstimatePx0==1)
                SES = diag(SEPx0terms);
            end
            if(PPEM_Constraints.Estimatex0==1)
                SEx0=SEx0terms;
            end

            if(PPEM_Constraints.EstimateA==1)
                if(PPEM_Constraints.AhatDiag==1)
                    SEA=diag(SEAterms);
                else
                    SEA=reshape(SEAterms,size(Ahat,2),size(Ahat,1))';
                end
            end
          
            if(PPEM_Constraints.QhatDiag==1)
                SEQ=diag(SEQterms);
            else
                SEQ=reshape(SEQterms,size(Qhat,2),size(Qhat,1))'; 
            end
            if(PPEM_Constraints.EstimateA==1)
                SE.A = SEA;
            end
            SE.Q = SEQ;
          
            if(PPEM_Constraints.EstimatePx0==1)
                SE.Px0=SES;
            end
            if(PPEM_Constraints.Estimatex0==1)
                SE.x0=SEx0;
            end
            
            SEMu = SEMuTerms;
            SEBeta=reshape(SEBetaTerms,size(betahat,2),size(betahat,1))';

            SE.mu = SEMu;
            SE.beta = SEBeta;
            if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                SEGamma=reshape(SEGammaTerms,size(gammahat,2),size(gammahat,1))';
                SE.gamma = SEGamma;
            end
            % Compute parameter p-values
            
            if(PPEM_Constraints.EstimateA==1)
                clear h p;
                if(PPEM_Constraints.AhatDiag==1)
                    VecParams = diag(Ahat);
                    VecSE     = diag(SEA);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pA = diag(p);
                else
                    VecParams = reshape(Ahat,[numel(Ahat) 1]);
                    VecSE     = reshape(SEA, [numel(Ahat) 1]);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end  
                    pA = reshape(p, [size(Ahat,1) size(Ahat,2)]);
                end
            end

            %Q matrix
            clear h p;
            if(PPEM_Constraints.QhatDiag==1)
                if(PPEM_Constraints.QhatIsotropic==1)
                    VecParams = Qhat(1,1);
                    VecSE     = SEQ(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pQ = diag(p);
                else
                    VecParams = diag(Qhat);
                    VecSE     = diag(SEQ);
                    for i=1:length(VecParams)
                       [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pQ = diag(p);
                end
            else
                VecParams = reshape(Qhat,[numel(Qhat) 1]);
                VecSE     = reshape(SEQ, [numel(Qhat) 1]);
                for i=1:length(VecParams)
                   [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end  
                pQ = reshape(p, [size(Qhat,1) size(Qhat,2)]);
            end
            %Px0
            if(PPEM_Constraints.EstimatePx0==1)
                clear h p;
                if(PPEM_Constraints.Px0Isotropic==1)
                    VecParams = Px0hat(1,1);
                    VecSE     = SES(1,1);
                    [h p] = ztest(VecParams,0,VecSE);
                    pPX0 = diag(p);
                else
                    VecParams = diag(Px0hat);
                    VecSE     = diag(SES);
                    for i=1:length(VecParams)
                        [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                    end
                    pPX0 = diag(p);
                end
            end


            if(PPEM_Constraints.Estimatex0==1)
                clear h p;
                VecParams = x0hat;
                VecSE     = SEx0;
                for i=1:length(VecParams)
                    [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end
                pX0 = p';
            end
            
            %Mu
            clear h p;
            VecParams = muhat;
            VecSE     = SEMu;
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end
            pMu = p';
            
            %Beta
            clear h p;
            VecParams = reshape(betahat,[numel(betahat),1]);
            VecSE     = reshape(SEBeta, [numel(SEBeta),1]);
            for i=1:length(VecParams)
                [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
            end    
            pBeta = reshape(p, [size(betahat,1) size(betahat,2)]);
            
            %Gamma
            clear h p;
            if((numel(gammahat)==1 && gammahat~=0) || numel(gammahat)>1)
                VecParams = reshape(gammahat,[numel(gammahat),1]);
                VecSE     = reshape(SEGamma, [numel(gammahat),1]);
                for i=1:length(VecParams)
                    [h(i) p(i)] = ztest(VecParams(i),0,VecSE(i));
                end    
                pGamma = reshape(p, [size(gammahat,1) size(gammahat,2)]);
            end
            if(PPEM_Constraints.EstimateA==1)
                Pvals.A = pA;
            end
            Pvals.Q = pQ;
           
            if(PPEM_Constraints.EstimatePx0==1)
                Pvals.Px0 = pPX0;
            end
            if(PPEM_Constraints.Estimatex0==1)
                Pvals.x0 = pX0;
            end
            Pvals.mu = pMu;
            Pvals.beta = pBeta;
            
            if(numel(gammahat)==1)
                if(gammahat~=0)
                    Pvals.gamma = pGamma;
                end
            else
                Pvals.gamma = pGamma;
            end

        end
        function [xKFinal,WKFinal,Ahat, Qhat, muhat, betahat, gammahat, x0hat, Px0hat, IC, SE, Pvals,nIter]=PP_EM(dN, Ahat0, Qhat0, mu, beta, fitType,delta, gamma, windowTimes, x0, Px0,PPEM_Constraints,MstepMethod)
            numStates = size(Ahat0,1);
            if(nargin<13 || isempty(MstepMethod))
               MstepMethod='GLM'; %or NewtonRaphson 
            end
            if(nargin<12 || isempty(PPEM_Constraints))
                PPEM_Constraints = DecodingAlgorithms.PP_EMCreateConstraints;
            end
            if(nargin<11 || isempty(Px0))
                Px0=10e-10*eye(numStates,numStates);
            end
            if(nargin<10 || isempty(x0))
                x0=zeros(numStates,1);
            end
            
            if(nargin<9 || isempty(windowTimes))
                if(isempty(gamma))
                    windowTimes =[];
                else
    %                 numWindows =length(gamma0)+1; 
                    windowTimes = 0:delta:(length(gamma)+1)*delta;
                end
            end
            if(nargin<8)
                gamma=[];
            end
            if(nargin<7 || isempty(delta))
                delta = .001;
            end
            if(nargin<6)
                fitType = 'poisson';
            end
            
            minTime=0;
            maxTime=(size(dN,2)-1)*delta;
            K=size(dN,1);
            if(~isempty(windowTimes))
                histObj = History(windowTimes,minTime,maxTime);
                for k=1:K
                    nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
                    nst{k}.setMinTime(minTime);
                    nst{k}.setMaxTime(maxTime);
%                     HkAll{k} = histObj.computeHistory(nst{k}).dataToMatrix;
                    HkAll(:,:,k) = histObj.computeHistory(nst{k}).dataToMatrix;
                end
            else
                for k=1:K
%                     HkAll{k} = 0;
                    HkAll(:,:,k) = 0;
                end
                gamma=0;
            end



    %         tol = 1e-3; %absolute change;
            tolAbs = 1e-3;
            tolRel = 1e-3;
            llTol  = 1e-3;
            cnt=1;

            maxIter = 100;

            
            A0 = Ahat0;
            Q0 = Qhat0;
           
            Ahat{1} = A0;
            Qhat{1} = Q0;
            x0hat{1} = x0;
            Px0hat{1} = Px0;
            muhat{1} = mu;
            betahat{1} = beta;
            gammahat{1} = gamma;
            numToKeep=10;
            scaledSystem=1;
            
            if(scaledSystem==1)
                Tq = eye(size(Qhat{1}))/(chol(Qhat{1}));
                Ahat{1}= Tq*Ahat{1}/Tq;
                Qhat{1}= Tq*Qhat{1}*Tq';
                x0hat{1} = Tq*x0;
                Px0hat{1} = Tq*Px0*Tq';
                betahat{1}=(betahat{1}'/Tq)';
            end

            cnt=1;
            dLikelihood(1)=inf;
%             x0hat = x0;
            negLL=0;
            IkedaAcc=PPEM_Constraints.EnableIkeda;
            %Forward EM
            stoppingCriteria =0;
%             logllNew= -inf;

            disp('                        Point-Process Observation EM Algorithm                        ');    
            while(stoppingCriteria~=1 && cnt<=maxIter)
                 storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
                 storeIndP1= mod(cnt,numToKeep)+1;
                 storeIndM1= mod(cnt-2,numToKeep)+1;
                disp('--------------------------------------------------------------------------------------------------------');
                disp(['Iteration #' num2str(cnt)]);
                disp('--------------------------------------------------------------------------------------------------------');
                

                [x_K{storeInd},W_K{storeInd},ll(cnt),ExpectationSums{storeInd}]=...
                    DecodingAlgorithms.PP_EStep(Ahat{storeInd},Qhat{storeInd},dN, muhat{storeInd}, betahat{storeInd},fitType,gammahat{storeInd},HkAll, x0hat{storeInd}, Px0hat{storeInd});
                
                [Ahat{storeIndP1}, Qhat{storeIndP1}, muhat{storeIndP1}, betahat{storeIndP1}, gammahat{storeIndP1},x0hat{storeIndP1},Px0hat{storeIndP1}] ...
                    = DecodingAlgorithms.PP_MStep(dN,x_K{storeInd},W_K{storeInd},x0hat{storeInd}, Px0hat{storeInd}, ExpectationSums{storeInd}, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,PPEM_Constraints,MstepMethod);
              
                if(IkedaAcc==1)
                    disp(['****Ikeda Acceleration Step****']);
                     
                     if(gammahat{storeIndP1}==0)% No history effect
                        dataMat = [ones(size(dN,2),1) x_K{storeInd}']; % design matrix: X 
                        coeffsMat = [muhat{storeIndP1} betahat{storeIndP1}']; % coefficient vector: beta
                        minTime=0;
                        maxTime=(size(dN,2)-1)*delta;
                        time=minTime:delta:maxTime;
                        clear nstNew;
                        for cc=1:length(muhat{storeIndP1})
                             tempData  = exp(dataMat*coeffsMat(cc,:)');

                             if(strcmp(fitType,'poisson'))
                                 lambdaData = tempData;
                             else
                                lambdaData = tempData./(1+tempData); % Conditional Intensity Function for ith cell
                             end
                             lambda{cc}=Covariate(time,lambdaData./delta, ...
                                 '\Lambda(t)','time','s','spikes/sec',...
                                 {strcat('\lambda_{',num2str(cc),'}')},{{' ''b'' '}});
                             lambda{cc}=lambda{cc}.resample(1/delta);

                             % generate one realization for each cell
                             tempSpikeColl{cc} = CIF.simulateCIFByThinningFromLambda(lambda{cc},1);          
                             nstNew{cc} = tempSpikeColl{cc}.getNST(1);     % grab the realization
                             nstNew{cc}.setName(num2str(cc));              % give each cell a unique name
%                              subplot(4,3,[8 11]);
%                              h2=lambda{cc}.plot([],{{' ''k'', ''LineWidth'' ,.5'}}); 
%                              legend off; hold all; % Plot the CIF

                        end
                        
                        spikeColl = nstColl(nstNew); % Create a neural spike train collection
                     else
                         time;
                     end
                     
                     dNNew=spikeColl.dataToMatrix';
                     dNNew(dNNew>1)=1; % more than one spike per bin will be treated as one spike. In
                                    % general we should pick delta small enough so that there is
                                    % only one spike per bin
                                    
                                    
%                                     [x_K,W_K,logll,ExpectationSums]=PP_EStep(A,Q,dN, mu, beta,fitType,gamma,HkAll, x0, Px0)
                     [x_KNew,W_KNew,logllNew,ExpectationSumsNew]=...
                        DecodingAlgorithms.PP_EStep(Ahat{storeInd},Qhat{storeInd},dNNew, muhat{storeInd}, betahat{storeInd},fitType,gammahat{storeInd},HkAll, x0, Px0);


                     [AhatNew, QhatNew, muhatNew, betahatNew, gammahatNew,x0new,Px0new] ...
                        = DecodingAlgorithms.PP_MStep(dNNew,x_KNew,W_KNew, x0hat{storeInd}, Px0hat{storeInd}, ExpectationSumsNew, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,PPEM_Constraints,MstepMethod);
               
                    Ahat{storeIndP1} = 2*Ahat{storeIndP1}-AhatNew;
                    Qhat{storeIndP1} = 2*Qhat{storeIndP1}-QhatNew;
                    Qhat{storeIndP1} = (Qhat{storeIndP1}+Qhat{storeIndP1}')/2;
                    muhat{storeIndP1}= 2*muhat{storeIndP1}-muhatNew;
                    betahat{storeIndP1} = 2*betahat{storeIndP1}-betahatNew;
                    gammahat{storeIndP1}= 2*gammahat{storeIndP1}-gammahatNew;
%                     x0hat{storeIndP1}   = 2*x0hat{storeIndP1} - x0new;
%                     Px0hat{storeIndP1}  = 2*Px0hat{storeIndP1}- Px0new;
%                     [V,D] = eig(Px0hat{storeIndP1});
%                     D(D<0)=1e-9;
%                     Px0hat{storeIndP1} = V*D*V';
%                     Px0hat{storeIndP1}  = (Px0hat{storeIndP1}+Px0hat{storeIndP1}')/2;
                    
               
                end
                if(PPEM_Constraints.EstimateA==0)
                    Ahat{storeIndP1}=Ahat{storeInd};
                end
                if(cnt==1)
                    dLikelihood(cnt+1)=inf;
                else
                    dLikelihood(cnt+1)=(ll(cnt)-ll(cnt-1));%./abs(ll(cnt-1));
                end
                if(cnt==1)
                    QhatInit = Qhat{1};
                    xKInit = x_K{1};
                end
                %Plot the progress
%                 if(mod(cnt,2)==0)
                if(cnt==1)
                    scrsz = get(0,'ScreenSize');
                    h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
                end
                figure(h);
                time = linspace(minTime,maxTime,size(x_K{storeInd},2));
                subplot(2,4,[1 2 5 6]); plot(1:cnt,ll,'k','Linewidth', 2); hy=ylabel('Log Likelihood'); hx=xlabel('Iteration'); axis auto;
                set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                subplot(2,4,3:4); hNew=plot(time, x_K{storeInd}','Linewidth', 2); hy=ylabel('States'); hx=xlabel('time [s]');
                set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold'); 
                hold on; hOrig=plot(time, xKInit','--','Linewidth', 2); 
                legend([hOrig(1) hNew(1)],'Initial','Current');
                  
                    
                subplot(2,4,7:8); hNew=plot(diag(Qhat{storeInd}),'o','Linewidth', 2); hy=ylabel('Q'); hx=xlabel('Diagonal Entry');
                set(gca, 'XTick'       , 1:1:length(diag(Qhat{storeInd})));
                set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
                hold on; hOrig=plot(diag(QhatInit),'r.','Linewidth', 2);
                legend([hOrig(1) hNew(1)],'Initial','Current');
                drawnow;
                hold off;

                if(cnt==1)
                    dMax=inf;
                else
                 dQvals = max(max(abs(sqrt(Qhat{storeInd})-sqrt(Qhat{storeIndM1}))));
                 dAvals = max(max(abs((Ahat{storeInd})-(Ahat{storeIndM1}))));
                 dMuvals = max(abs((muhat{storeInd})-(muhat{storeIndM1})));
                 dBetavals = max(max(abs((betahat{storeInd})-(betahat{storeIndM1}))));
                 dGammavals = max(max(abs((gammahat{storeInd})-(gammahat{storeIndM1}))));
                 dMax = max([dQvals,dAvals,dMuvals,dBetavals,dGammavals]);
                end

                if(cnt==1)
                    disp(['Max Parameter Change: N/A']);
                else
                    disp(['Max Parameter Change: ' num2str(dMax)]);
                end
                
                cnt=(cnt+1);
                if(dMax<tolAbs)
                    stoppingCriteria=1;
                    display(['         EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
                    negLL=0;
                end
            
                if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
                    stoppingCriteria=1;
                    display(['         EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
                    negLL=1;
                end
                

            end
            disp('--------------------------------------------------------------------------------------------------------');

            maxLLIndex  = find(ll == max(ll),1,'first');
            maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
            if(maxLLIndex==1)
%                 maxLLIndex=cnt-1;
                maxLLIndex =1;
                maxLLIndMod = 1;
            elseif(isempty(maxLLIndex))
               maxLLIndex = 1; 
               maxLLIndMod = 1;
%             else
%                maxLLIndMod = mod(maxLLIndex,numToKeep); 
               
            end
            nIter   = cnt-1;  
%             maxLLIndMod
           
            xKFinal = x_K{maxLLIndMod};
            WKFinal = W_K{maxLLIndMod};
            Ahat = Ahat{maxLLIndMod};
            Qhat = Qhat{maxLLIndMod};
            muhat= muhat{maxLLIndMod};
            betahat = betahat{maxLLIndMod};
            gammahat = gammahat{maxLLIndMod};
            x0hat =x0hat{maxLLIndMod};
            Px0hat=Px0hat{maxLLIndMod};
            
             if(scaledSystem==1)
               Tq = eye(size(Qhat))/(chol(Q0));
               Ahat=Tq\Ahat*Tq;
               Qhat=(Tq\Qhat)/Tq';
               xKFinal = Tq\xKFinal;
               x0hat = Tq\x0hat;
               Px0hat= (Tq\Px0hat)/(Tq');
               tempWK =zeros(size(WKFinal));
               for kk=1:size(WKFinal,3)
                tempWK(:,:,kk)=(Tq\WKFinal(:,:,kk))/Tq';
               end
               WKFinal = tempWK;
               betahat=(betahat'*Tq)';
             end
            llFinal=ll(end);
            ll = ll(maxLLIndex);
            ExpectationSumsFinal = ExpectationSums{maxLLIndMod};
            if(nargout>10)
                [SE, Pvals]=DecodingAlgorithms.PP_ComputeParamStandardErrors(dN,...
                    xKFinal, WKFinal, Ahat, Qhat, x0hat, Px0hat, ExpectationSumsFinal,...
                    fitType, muhat, betahat, gammahat, windowTimes, HkAll,...
                    PPEM_Constraints);
            end

             %Compute number of parameters
            if(PPEM_Constraints.EstimateA==1 && PPEM_Constraints.AhatDiag==1)
                n1=size(Ahat,1); 
            elseif(PPEM_Constraints.EstimateA==1 && PPEM_Constraints.AhatDiag==0)
                n1=numel(Ahat);
            else 
                n1=0;
            end
            if(PPEM_Constraints.QhatDiag==1 && PPEM_Constraints.QhatIsotropic==1)
                n2=1;
            elseif(PPEM_Constraints.QhatDiag==1 && PPEM_Constraints.QhatIsotropic==0)
                n2=size(Qhat,1);
            else
                n2=numel(Qhat);
            end


            if(PPEM_Constraints.EstimatePx0==1 && PPEM_Constraints.Px0Isotropic==1)
                n3=1;
            elseif(PPEM_Constraints.EstimatePx0==1 && PPEM_Constraints.Px0Isotropic==0)
                n3=size(Px0hat,1);
            else
                n3=0;
            end

            if(PPEM_Constraints.Estimatex0==1)   
                n4=size(x0hat,1);
            else
                n4=0;
            end

            n5=size(muhat,1);
            n6=numel(betahat);
            if(numel(gammahat)==1)
                if(gammahat==0)
                    n7=0;
                else
                    n7=1;
                end
            else
                n7=numel(gammahat);
            end
            nTerms=n1+n2+n3+n4+n5+n6+n7;
            
            K  = size(xKFinal,2); 
            Dx = size(Ahat,2);
            sumXkTerms = ExpectationSums{maxLLIndMod}.sumXkTerms;
            llobs = ll + Dx*K/2*log(2*pi)+K/2*log(det(Qhat))...
                + 1/2*trace(Qhat\sumXkTerms)...
                + Dx/2*log(2*pi)+1/2*log(det(Px0hat)) ...
                + 1/2*Dx;
            AIC = 2*nTerms - 2*llobs;
            AICc= AIC+ 2*nTerms*(nTerms+1)/(K-nTerms-1);
            BIC = -2*llobs+nTerms*log(K);
            IC.AIC = AIC;
            IC.AICc= AICc;
            IC.BIC = BIC;
            IC.llobs = llobs;
            IC.llcomp=ll;
           
            
        end
        %         function  [xKFinal,WKFinal,Ahat, Qhat, muhat, betahat, gammahat, x0hat, Px0hat, logll,nIter,negLL]=PP_EM(dN, Ahat0, Qhat0, mu, beta, fitType,delta, gamma, windowTimes, x0, Px0,MstepMethod)
%             numStates = size(Ahat0,1);
%             if(nargin<12 || isempty(MstepMethod))
%                MstepMethod='GLM'; %or NewtonRaphson 
%             end
%             if(nargin<11 || isempty(Px0))
%                 Px0=10e-10*eye(numStates,numStates);
%             end
%             if(nargin<10 || isempty(x0))
%                 x0=zeros(numStates,1);
%             end
%             
%             if(nargin<9 || isempty(windowTimes))
%                 if(isempty(gamma)||gamma==0)
%                     windowTimes =[];
%                 else
%     %                 numWindows =length(gamma0)+1; 
%                     windowTimes = 0:delta:(length(gamma)+1)*delta;
%                 end
%             end
%             if(nargin<8)
%                 gamma=[];
%             end
%             if(nargin<11 || isempty(delta))
%                 delta = .001;
%             end
%             if(nargin<6)
%                 fitType = 'poisson';
%             end
%             
%             minTime=0;
%             maxTime=(size(dN,2)-1)*delta;
%             K=size(dN,1);
%             N=size(dN,2);
%             if(~isempty(windowTimes))
%                 histObj = History(windowTimes,minTime,maxTime);
%                 for k=1:K
%                     nst{k} = nspikeTrain( (find(dN(k,:)==1)-1)*delta);
%                     nst{k}.setMinTime(minTime);
%                     nst{k}.setMaxTime(maxTime);
% %                     HkAll{k} = histObj.computeHistory(nst{k}).dataToMatrix;
%                     HkAll(:,:,k) = histObj.computeHistory(nst{k}).dataToMatrix;
%                 end
%                 if(size(gamma,1)==K)
%                     gamma=gamma';
%                 end
%                 
%             else
%                 for k=1:K
%                     HkAll(:,:,k) = zeros(N,length(windowTimes)-1);
%                 end
%                 gamma=0;
%             end
%                 
% 
% 
%     %         tol = 1e-3; %absolute change;
%             tolAbs = 1e-3;
%             tolRel = 1e-3;
%             llTol  = 1e-3;
%             cnt=1;
% 
%             maxIter = 100;
% 
%             
%             A0 = Ahat0;
%             Q0 = Qhat0;
%            
%             Ahat{1} = A0;
%             Qhat{1} = Q0;
%             x0hat{1} = x0;
%             Px0hat{1} = Px0;
%             muhat{1} = mu;
%             betahat{1} = beta;
%             gammahat{1} = gamma;
%             numToKeep=10;
%             scaledSystem=1;
%             
%             if(scaledSystem==1)
%                 Tq = eye(size(Qhat{1}))/(chol(Qhat{1}));
%                 Ahat{1}= Tq*Ahat{1}/Tq;
%                 Qhat{1}= Tq*Qhat{1}*Tq';
%                 x0hat{1} = Tq*x0;
%                 Px0hat{1} = Tq*Px0*Tq';
%                 betahat{1}=(betahat{1}'/Tq)';
%             end
% 
%             cnt=1;
%             dLikelihood(1)=inf;
%             negLL=0;
%             IkedaAcc=0;
%             %Forward EM
%             stoppingCriteria =0;
%                 
%             while(stoppingCriteria~=1 && cnt<=maxIter)
%                  storeInd = mod(cnt-1,numToKeep)+1; %make zero-based then mod, then add 1
%                  storeIndP1= mod(cnt,numToKeep)+1;
%                  storeIndM1= mod(cnt-2,numToKeep)+1;
%                 disp('---------------');
%                 disp(['Iteration #' num2str(cnt)]);
%                 disp('---------------');
%                 
%                 
%                 [x_K{storeInd},W_K{storeInd},logll(cnt),ExpectationSums{storeInd}]=...
%                     DecodingAlgorithms.PP_EStep(Ahat{storeInd},Qhat{storeInd},dN, muhat{storeInd}, betahat{storeInd},fitType,gammahat{storeInd},HkAll, x0hat{storeInd}, Px0hat{storeInd});
%                 
%                 [Ahat{storeIndP1}, Qhat{storeIndP1}, muhat{storeIndP1}, betahat{storeIndP1}, gammahat{storeIndP1},x0hat{storeIndP1},Px0hat{storeIndP1}] ...
%                     = DecodingAlgorithms.PP_MStep(dN,x_K{storeInd},W_K{storeInd},x0hat{storeInd},ExpectationSums{storeInd}, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,MstepMethod);
%               
%                 if(IkedaAcc==1)
%                     disp(['****Ikeda Acceleration Step****']);
%                      
%                      if(gammahat{storeIndP1}==0)% No history effect
%                         dataMat = [ones(size(dN,2),1) x_K{storeInd}']; % design matrix: X 
%                         coeffsMat = [muhat{storeIndP1} betahat{storeIndP1}']; % coefficient vector: beta
%                         minTime=0;
%                         maxTime=(size(dN,2)-1)*delta;
%                         time=minTime:delta:maxTime;
%                         clear nstNew;
%                         for cc=1:length(muhat{storeIndP1})
%                              tempData  = exp(dataMat*coeffsMat(cc,:)');
% 
%                              if(strcmp(fitType,'poisson'))
%                                  lambdaData = tempData;
%                              else
%                                 lambdaData = tempData./(1+tempData); % Conditional Intensity Function for ith cell
%                              end
%                              lambda{cc}=Covariate(time,lambdaData./delta, ...
%                                  '\Lambda(t)','time','s','spikes/sec',...
%                                  {strcat('\lambda_{',num2str(cc),'}')},{{' ''b'' '}});
%                              lambda{cc}=lambda{cc}.resample(1/delta);
% 
%                              % generate one realization for each cell
%                              tempSpikeColl{cc} = CIF.simulateCIFByThinningFromLambda(lambda{cc},1);          
%                              nstNew{cc} = tempSpikeColl{cc}.getNST(1);     % grab the realization
%                              nstNew{cc}.setName(num2str(cc));              % give each cell a unique name
% %                              subplot(4,3,[8 11]);
% %                              h2=lambda{cc}.plot([],{{' ''k'', ''LineWidth'' ,.5'}}); 
% %                              legend off; hold all; % Plot the CIF
% 
%                         end
%                         
%                         spikeColl = nstColl(nstNew); % Create a neural spike train collection
%                      else
%                          time;
%                      end
%                      
%                      dNNew=spikeColl.dataToMatrix';
%                      dNNew(dNNew>1)=1; % more than one spike per bin will be treated as one spike. In
%                                     % general we should pick delta small enough so that there is
%                                     % only one spike per bin
%                                     
%                                     
% %                                     [x_K,W_K,logll,ExpectationSums]=PP_EStep(A,Q,dN, mu, beta,fitType,gamma,HkAll, x0, Px0)
%                      [x_KNew,W_KNew,logllNew,ExpectationSumsNew]=...
%                         DecodingAlgorithms.PP_EStep(Ahat{storeInd},Qhat{storeInd},dNNew, muhat{storeInd}, betahat{storeInd},fitType,gammahat{storeInd},HkAll, x0, Px0);
% 
%                 
%                      [AhatNew, QhatNew, muhatNew, betahatNew, gammahatNew,x0new,Px0new] ...
%                         = DecodingAlgorithms.PP_MStep(dNNew,x_KNew,W_KNew, x0hat{storeInd}, ExpectationSumsNew, fitType,muhat{storeInd},betahat{storeInd}, gammahat{storeInd},windowTimes,HkAll,MstepMethod);
%                
%                     Ahat{storeIndP1} = 2*Ahat{storeIndP1}-AhatNew;
%                     Qhat{storeIndP1} = 2*Qhat{storeIndP1}-QhatNew;
%                     Qhat{storeIndP1} = (Qhat{storeIndP1}+Qhat{storeIndP1}')/2;
%                     muhat{storeIndP1}= 2*muhat{storeIndP1}-muhatNew;
%                     betahat{storeIndP1} = 2*betahat{storeIndP1}-betahatNew;
%                     gammahat{storeIndP1}= 2*gammahat{storeIndP1}-gammahatNew;
% %                     x0hat{storeIndP1}   = 2*x0hat{storeIndP1} - x0new;
% %                     Px0hat{storeIndP1}  = 2*Px0hat{storeIndP1}- Px0new;
% %                     [V,D] = eig(Px0hat{storeIndP1});
% %                     D(D<0)=1e-9;
% %                     Px0hat{storeIndP1} = V*D*V';
% %                     Px0hat{storeIndP1}  = (Px0hat{storeIndP1}+Px0hat{storeIndP1}')/2;
%                     
%                
%                 end
% 
%                 if(cnt==1)
%                     dLikelihood(cnt+1)=inf;
%                 else
%                     dLikelihood(cnt+1)=(logll(cnt)-logll(cnt-1));%./abs(logll(cnt-1));
%                 end
%                 if(cnt==1)
%                     QhatInit = Qhat{1};
%                     xKInit = x_K{1};
%                 end
%                 %Plot the progress
% %                 if(mod(cnt,2)==0)
%                 if(cnt==1)
%                     scrsz = get(0,'ScreenSize');
%                     h=figure('OuterPosition',[scrsz(3)*.01 scrsz(4)*.04 scrsz(3)*.98 scrsz(4)*.95]);
%                 end
%                     figure(h);
%                     time = linspace(minTime,maxTime,size(x_K{storeInd},2));
%                     subplot(2,4,[1 2 5 6]); plot(1:cnt,logll,'k','Linewidth', 2); hy=ylabel('Log Likelihood'); hx=xlabel('Iteration'); axis auto;
%                     set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
%                     subplot(2,4,3:4); hNew=plot(time, x_K{storeInd}','Linewidth', 2); hy=ylabel('States'); hx=xlabel('time [s]');
%                     set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold'); 
%                     hold on; hOrig=plot(time, xKInit','--','Linewidth', 2); 
%                     legend([hOrig(1) hNew(1)],'Initial','Current');
%                   
%                     
%                     subplot(2,4,7:8); hNew=plot(diag(Qhat{storeInd}),'o','Linewidth', 2); hy=ylabel('Q'); hx=xlabel('Diagonal Entry');
%                     set(gca, 'XTick'       , 1:1:length(diag(Qhat{storeInd})));
%                     set([hx, hy],'FontName', 'Arial','FontSize',12,'FontWeight','bold');
%                     hold on; hOrig=plot(diag(QhatInit),'r.','Linewidth', 2);
%                     legend([hOrig(1) hNew(1)],'Initial','Current');
%                     drawnow;
%                     hold off;
% %                 end
%                 
%                 if(cnt==1)
%                     dMax=inf;
%                 else
%                  dQvals = max(max(abs(sqrt(Qhat{storeInd})-sqrt(Qhat{storeIndM1}))));
%                  dAvals = max(max(abs((Ahat{storeInd})-(Ahat{storeIndM1}))));
%                  dMuvals = max(abs((muhat{storeInd})-(muhat{storeIndM1})));
%                  dBetavals = max(max(abs((betahat{storeInd})-(betahat{storeIndM1}))));
%                  dGammavals = max(max(abs((gammahat{storeInd})-(gammahat{storeIndM1}))));
%                  dMax = max([dQvals,dAvals,dMuvals,dBetavals,dGammavals]);
%                 end
% 
% % 
% %                 dQRel = max(abs(dQvals./sqrt(Qhat(:,storeIndM1))));
% %                 dGammaRel = max(abs(dGamma./gammahat(storeIndM1,:)));
% %                 dMaxRel = max([dQRel,dGammaRel]);
% 
%                  
%                 cnt=(cnt+1);
%                 if(dMax<tolAbs)
%                     stoppingCriteria=1;
%                     display(['         EM converged at iteration# ' num2str(cnt-1) ' b/c change in params was within criteria']);
%                     negLL=0;
%                 end
%             
%                 if(abs(dLikelihood(cnt))<llTol  || dLikelihood(cnt)<0)
%                     stoppingCriteria=1;
%                     display(['         EM stopped at iteration# ' num2str(cnt-1) ' b/c change in likelihood was negative']);
%                     negLL=1;
%                 end
%                 
% 
%             end
%             
%             
% 
% 
%             maxLLIndex  = find(logll == max(logll),1,'first');
%             maxLLIndMod =  mod(maxLLIndex-1,numToKeep)+1;
%             if(maxLLIndex==1)
% %                 maxLLIndex=cnt-1;
%                 maxLLIndex =1;
%                 maxLLIndMod = 1;
%             elseif(isempty(maxLLIndex))
%                maxLLIndex = 1; 
%                maxLLIndMod = 1;
% %             else
% %                maxLLIndMod = mod(maxLLIndex,numToKeep); 
%                
%             end
%             nIter   = cnt-1;  
% %             maxLLIndMod
%            
%             xKFinal = x_K{maxLLIndMod};
%             WKFinal = W_K{maxLLIndMod};
%             Ahat = Ahat{maxLLIndMod};
%             Qhat = Qhat{maxLLIndMod};
%             muhat= muhat{maxLLIndMod};
%             betahat = betahat{maxLLIndMod};
%             gammahat = gammahat{maxLLIndMod};
%             x0hat =x0hat{maxLLIndMod};
%             Px0hat=Px0hat{maxLLIndMod};
%             
%              if(scaledSystem==1)
%                Tq = eye(size(Qhat))/(chol(Q0));
%                Ahat=Tq\Ahat*Tq;
%                Qhat=(Tq\Qhat)/Tq';
%                xKFinal = Tq\xKFinal;
%                x0hat = Tq\x0hat;
%                Px0hat= (Tq\Px0hat)/(Tq');
%                tempWK =zeros(size(WKFinal));
%                for kk=1:size(WKFinal,3)
%                 tempWK(:,:,kk)=(Tq\WKFinal(:,:,kk))/Tq';
%                end
%                WKFinal = tempWK;
%                betahat=(betahat'*Tq)';
%              end
%             
%             logll = logll(maxLLIndex);
%             ExpectationSumsFinal = ExpectationSums{maxLLIndMod};
%             K=size(dN,1);
%             SumXkTermsFinal = diag(Qhat(:,:,end))*K;
%             logllFinal=logll(end);
%             McInfo=100;
%             McCI = 3000;
% 
% %             nIter = [];%[nIter1,nIter2,nIter3];
%   
%             
%             K  = size(dN,1); 
%             Dx = size(Ahat,2);
%             sumXkTerms = ExpectationSums{maxLLIndMod}.sumXkTerms;
%             logllobs = logll + Dx*K/2*log(2*pi)+K/2*log(det(Qhat))+ 1/2*trace(pinv(Qhat)*sumXkTerms); 
%                   
% %             InfoMat = DecodingAlgorithms.estimateInfoMat_mPPCO(fitType,xKFinal, WKFinal,Ahat,Qhat,Chat, Rhat,alphahat, muhat, betahat,gammahat,dN,windowTimes, HkAll,delta,ExpectationSums{maxLLIndMod},McInfo);
% %             
% %             
% %             fitResults = DecodingAlgorithms.prepareEMResults(fitType,neuronName,dN,HkAll,xKFinal,WKFinal,Qhat,gammahat,windowTimes,delta,InfoMat,logllobs);
% %             [stimCIs, stimulus] = DecodingAlgorithms.ComputeStimulusCIs(fitType,xKFinal,WkuFinal,delta,McCI);
% %             
%            
%         end
        function [x_K,W_K,logll,ExpectationSums]=PP_EStep(A,Q,dN, mu, beta,fitType,gamma,HkAll, x0, Px0)
            DEBUG = 0;
            [numCells,K]   = size(dN); 
            Dx = size(A,2);
            
            x_p     = zeros( size(A,2), K+1 );
            x_u     = zeros( size(A,2), K );
            W_p    = zeros( size(A,2),size(A,2), K+1 );
            W_u    = zeros( size(A,2),size(A,2), K );
            x_p(:,1)= A(:,:)*x0;
            W_p(:,:,1)=A*Px0*A' + Q;
            HkPerm=permute(HkAll, [2 3 1]);
            for k=1:K
                [x_u(:,k), W_u(:,:,k)] = DecodingAlgorithms.PPDecode_updateLinear(x_p(:,k), W_p(:,:,k), dN,mu,beta,fitType,gamma,HkPerm,k,[]);
                [x_p(:,k+1), W_p(:,:,k+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,k), W_u(:,:,k), A(:,:,min(size(A,3),k)), Q(:,:,min(size(Q,3))));
            end  
            
            [x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u); 
             
            numStates = size(x_K,1);
            Wku=zeros(numStates,numStates,K,K);
            Tk = zeros(numStates,numStates,K-1);
            for k=1:K
                Wku(:,:,k,k)=W_K(:,:,k);
            end

            for u=K:-1:2
                for k=(u-1):-1:(u-1)
                    Tk(:,:,k)=A;
%                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
                    Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
                    Wku(:,:,u,k)=Wku(:,:,k,u)';
                end
            end
            
            %All terms
            Sxkm1xk = zeros(Dx,Dx);
            Sxkm1xkm1 = zeros(Dx,Dx);
            Sxkxk = zeros(Dx,Dx);
            for k=1:K
                if(k==1)
                    Sxkm1xk   = Sxkm1xk+Px0*A'/W_p(:,:,1)*Wku(:,:,1,1);
                    Sxkm1xkm1 = Sxkm1xkm1+Px0+x0*x0';     
                else
                      Sxkm1xk =  Sxkm1xk+Wku(:,:,k-1,k)+x_K(:,k-1)*x_K(:,k)';
                      Sxkm1xkm1= Sxkm1xkm1+Wku(:,:,k-1,k-1)+x_K(:,k-1)*x_K(:,k-1)';
                end
                Sxkxk = Sxkxk+Wku(:,:,k,k)+x_K(:,k)*x_K(:,k)';

            end
            Sxkxk = 0.5*(Sxkxk+Sxkxk');
            sumXkTerms = Sxkxk-A*Sxkm1xk-Sxkm1xk'*A'+A*Sxkm1xkm1*A';
            Sxkxkm1 = Sxkm1xk';
            
            %Vectorize for loop over cells
            if(strcmp(fitType,'poisson'))
                sumPPll=0;
                Histtermperm = permute(HkAll,[2 3 1]);
                
                for k=1:K
%                    Hk=squeeze(HkAll(k,:,:)); 
                   Hk= Histtermperm(:,:,k);
                   if(size(Hk,1)==numCells)
                       Hk = Hk';
                   end
                   xk = x_K(:,k);
                   if(numel(gamma)==1)
                        gammaC=repmat(gamma,1,numCells);
                   else 
                        gammaC=gamma;
                   end
                   terms=mu+beta'*xk+diag(gammaC'*Hk);
                   Wk = W_K(:,:,k);
                   ld = exp(terms);
                   bt = beta;
                   ExplambdaDelta =ld+0.5*(ld.*diag((bt'*Wk*bt)));
                   ExplogLD = terms;
                   sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta);
                        
                end
                
            %Vectorize over number of cells
            elseif(strcmp(fitType,'binomial'))
                sumPPll=0;
                Histtermperm = permute(HkAll,[2 3 1]);
                for k=1:K
%                     Hk=squeeze(HkAll(k,:,:)); 
                    Hk= Histtermperm(:,:,k);
                    if(size(Hk,1)==numCells)
                       Hk = Hk';
                    end
                    xk = x_K(:,k);
                    if(numel(gamma)==1)
                        gammaC=repmat(gamma,1,numCells);
                    else 
                        gammaC=gamma;
                   end
                   terms=mu+beta'*xk+diag(gammaC'*Hk);
                   Wk = W_K(:,:,k);
                   ld = exp(terms)./(1+exp(terms));
                   bt = beta;     
                   ExplambdaDelta = ld+0.5*(ld.*(1-ld).*(1-2.*ld)).*diag((bt'*Wk*bt));
                   ExplogLD = log(ld)+0.5*(-ld.*(1-ld)).*diag(bt'*Wk*bt);
                   sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta); 
                    
                end

                
            end

            logll = -Dx*K/2*log(2*pi)-K/2*log(det(Q)) ...
                    - Dx/2*log(2*pi) -1/2*log(det(Px0))  ...
                    +sumPPll - 1/2*trace((eye(size(Q))/Q)*sumXkTerms) ...
                    -Dx/2;
                string0 = ['logll: ' num2str(logll)];
                disp(string0);
                if(DEBUG==1)
                    string1 = ['-K/2*log(det(Q)):' num2str(-K/2*log(det(Q)))];
                    string12= ['Constants: ' num2str(-Dx*K/2*log(2*pi)- Dx/2*log(2*pi) -Dx/2 -1/2*log(det(Px0)))];
                    string2 = ['SumPPll: ' num2str(sumPPll)];
                    string3 = ['-.5*trace(Q\sumXkTerms): ' num2str(-.5*trace(Q\sumXkTerms))];
                   
                    disp(string1);
                    disp(['Q=' num2str(diag(Q)')]);
                    disp(string12);
                    disp(string2);
                    disp(string3);
                end

                ExpectationSums.Sxkm1xkm1=Sxkm1xkm1;
                ExpectationSums.Sxkm1xk=Sxkm1xk;
                ExpectationSums.Sxkxkm1=Sxkxkm1;
                ExpectationSums.Sxkxk=Sxkxk;
                ExpectationSums.sumXkTerms=sumXkTerms;
                ExpectationSums.sumPPll=sumPPll;

        end
        %         function [x_K,W_K,logll,ExpectationSums]=PP_EStep(A,Q,dN, mu, beta,fitType,gamma,HkAll, x0, Px0)
%              
%             DEBUG = 0;
%             [numCells,K]   = size(dN); 
%             Dx = size(A,2);
%             
%             x_p     = zeros( size(A,2), K+1 );
%             x_u     = zeros( size(A,2), K );
%             W_p    = zeros( size(A,2),size(A,2), K+1 );
%             W_u    = zeros( size(A,2),size(A,2), K );
%             x_p(:,1)= A(:,:)*x0;
%             W_p(:,:,1)=A*Px0*A' + Q;
% %             WuConv=[];
%             for k=1:K
%                 [x_u(:,k), W_u(:,:,k)] = DecodingAlgorithms.PPDecode_updateLinear(x_p(:,k), W_p(:,:,k), dN,mu,beta,fitType,gamma,HkAll,k,[]);
%                 [x_p(:,k+1), W_p(:,:,k+1)] = DecodingAlgorithms.PPDecode_predict(x_u(:,k), W_u(:,:,k), A(:,:,min(size(A,3),k)), Q(:,:,min(size(Q,3))));
% %                 if(k>1 && isempty(WuConv))
% %                     diffWu = abs(W_u(:,:,k)-W_u(:,:,k-1));
% %                     maxWu  = max(max(diffWu));
% %                     if(maxWu<5e-2)
% %                         WuConv = W_u(:,:,k);
% %                         WuConvIter = k;
% %                     end
% %                 end
%             end
%      
%             
%             [x_K, W_K,Lk] = DecodingAlgorithms.kalman_smootherFromFiltered(A, x_p, W_p, x_u, W_u);
%             
%             %Best estimates of initial states given the data
%             W1G0 = A*Px0*A' + Q;
%             L0=Px0*A'/W1G0;
%             
%             Ex0Gy = x0+L0*(x_K(:,1)-x_p(:,1));        
%             Px0Gy = Px0+L0*(eye(size(W_K(:,:,1)))/(W_K(:,:,1))-eye(size(W1G0))/W1G0)*L0';
%             Px0Gy = (Px0Gy+Px0Gy')/2;
%             numStates = size(x_K,1);
%             Wku=zeros(numStates,numStates,K,K);
%             Tk = zeros(numStates,numStates,K-1);
%             for k=1:K
%                 Wku(:,:,k,k)=W_K(:,:,k);
%             end
% 
%             for u=K:-1:2
%                 for k=(u-1):-1:(u-1)
%                     Tk(:,:,k)=A;
% %                     Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'*pinv(W_p(:,:,k)); %From deJong and MacKinnon 1988
%                      Dk(:,:,k)=W_u(:,:,k)*Tk(:,:,k)'/(W_p(:,:,k+1)); %From deJong and MacKinnon 1988
%                     Wku(:,:,k,u)=Dk(:,:,k)*Wku(:,:,k+1,u);
%                     Wku(:,:,u,k)=Wku(:,:,k,u)';
%                 end
%             end
%             
%             %All terms
%             Sxkm1xk = zeros(Dx,Dx);
%             Sxkxkm1 = zeros(Dx,Dx);
%             Sxkm1xkm1 = zeros(Dx,Dx);
%             Sxkxk = zeros(Dx,Dx);
%          
%             for k=1:K
%                 if(k==1)
%                     Sxkm1xk   = Sxkm1xk+Px0*A'/W_p(:,:,1)*Wku(:,:,1,1);
%                     Sxkm1xkm1 = Sxkm1xkm1+Px0+x0*x0';     
%                 else
% %                   
%                       Sxkm1xk =  Sxkm1xk+Wku(:,:,k-1,k)+x_K(:,k-1)*x_K(:,k)';
%                        
%                       Sxkm1xkm1= Sxkm1xkm1+Wku(:,:,k-1,k-1)+x_K(:,k-1)*x_K(:,k-1)';
%                 end
%                 Sxkxk = Sxkxk+Wku(:,:,k,k)+x_K(:,k)*x_K(:,k)';
%                 
%             end
%             Sx0x0 = Px0+x0*x0';
%             Sxkxk = 0.5*(Sxkxk+Sxkxk');
%             sumXkTerms = Sxkxk-A*Sxkm1xk-Sxkm1xk'*A'+A*Sxkm1xkm1*A';
%             Sxkxkm1 = Sxkm1xk';
%             
% %             if(strcmp(fitType,'poisson'))
% %                 sumPPll=0;
% %                 for c=1:numCells
% % %                     Hk=HkAll{c};
% %                     Hk=squeeze(HkAll(k,:,c));
% %                     for k=1:K
% %                         xk = x_K(:,k);
% %                         if(numel(gamma)==1)
% %                             gammaC=gamma;
% %                         else 
% %                             gammaC=gamma(:,c);
% %                         end
% % %                         terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
% %                         if(numel(Hk)~=1)
% %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
% %                         else
% %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk';
% %                         end
% %                         Wk = W_K(:,:,k);
% %                         ld = exp(terms);
% %                         bt = beta(:,c);
% %                         ExplambdaDelta =ld+0.5*trace(bt*bt'*ld*Wk);
% %                         ExplogLD = terms;
% %                         sumPPll=sumPPll+dN(c,k).*ExplogLD - ExplambdaDelta;
% %                     end
% %                   
% %                             
% %                 end
% %             elseif(strcmp(fitType,'binomial'))
% %                 sumPPll=0;
% %                 for c=1:C
% %                     for k=1:K
% %                         Hk=squeeze(HkAll(k,:,c));
% %                         xk = x_K(:,k);
% %                         if(numel(gamma)==1)
% %                             gammaC=gamma;
% %                         else 
% %                             gammaC=gamma(:,c);
% %                         end
% %                         if(numel(Hk)~=1)
% %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
% %                         else
% %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk';
% %                         end
% %                         Wk = W_K(:,:,k);
% %                         ld = exp(terms)./(1+exp(terms));
% %                         bt = beta;
% %                         ExplambdaDelta =sum(ld+0.5*sum(bt'*bt*repmat(ld.*(1-ld).*(1-2.*ld),1,2)*Wk,2));
% %                         ExplogLD = (log(ld)+0.5*sum(bt*bt'*(repmat(ld.*(1-ld),1,size(bt,1))*Wk)')');
% %                         sumPPll=sumPPll+dN(:,k)'*ExplogLD - ExplambdaDelta;
% %                     end
% %                 end
% %                 
% % %                 for c=1:numCells
% % %                     Hk=HkAll{c};
% % %                     for k=1:K
% % %                         xk = x_K(:,k);
% % %                         if(numel(gamma)==1)
% % %                             gammaC=gamma;
% % %                         else 
% % %                             gammaC=gamma(:,c);
% % %                         end
% % %                         if(numel(Hk)~=1)
% % %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk(k,:)';
% % %                         else
% % %                             terms=mu(c)+beta(:,c)'*xk+gammaC'*Hk';
% % %                         end
% % %                         Wk = W_K(:,:,k);
% % %                         ld = exp(terms)./(1+exp(terms));
% % %                         bt = beta(:,c);
% % %                         ExplambdaDelta =ld+0.5*trace(bt*bt'*ld*(1-ld)*(1-2*ld)*Wk);
% % %                         ExplogLD = log(ld)+0.5*trace(-(bt*bt'*ld*(1-ld))*Wk);
% % %                         sumPPll=sumPPll+dN(c,k).*ExplogLD - ExplambdaDelta;
% % %                     end
% % %                   
% % %                             
% % %                 end
% %             end
% 
%             %Vectorize for loop over cells
%             if(strcmp(fitType,'poisson'))
%                 sumPPll=0;
%                 for k=1:K
%                    Hk=squeeze(HkAll(k,:,:)); 
%                    if(size(Hk,1)==numCells)
%                     Hk = Hk';
%                    end
%                    xk = x_K(:,k);
%                    if(numel(gamma)==1)
%                         gammaC=repmat(gamma,1,numCells);
%                    else 
%                         gammaC=gamma;
%                    end
% %                    if(size(gammaC,1)~=size(mu,1))
% %                         gammaC = gammaC';
% %                     end
% %                     if(size(Hk,1)~=size(mu,1))
% %                         Hk=Hk';
% %                     end
%                    terms=mu+beta'*xk+diag(gammaC'*Hk);
%                    Wk = W_K(:,:,k);
%                    ld = exp(terms);
%                    bt = beta;
%                    ExplambdaDelta =ld+0.5*(ld.*diag((bt'*Wk*bt)));
%                    ExplogLD = terms;
%                    sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta);
%                         
%                 end
%                 
%             %Vectorize over number of cells
%             elseif(strcmp(fitType,'binomial'))
%                 sumPPll=0;
%                 for k=1:K
%                     Hk=squeeze(HkAll(k,:,:));
%                     if(size(Hk,1)==numCells)
%                        Hk = Hk';
%                     end
%                     xk = x_K(:,k);
%                     if(numel(gamma)==1)
%                         gammaC=repmat(gamma,1,numCells);
%                     else 
%                         gammaC=gamma;
%                     end
% %                     if(size(gammaC,1)~=size(mu,1))
% %                         gammaC = gammaC';
% %                     end
% %                     if(size(Hk,1)~=size(mu,1))
% %                         Hk=Hk';
% %                     end
%                    terms=mu+beta'*xk+diag(gammaC'*Hk);
%                    Wk = W_K(:,:,k);
%                    ld = exp(terms)./(1+exp(terms));
%                    bt = beta;     
%                    ExplambdaDelta = ld+0.5*(ld.*(1-ld).*(1-2.*ld)).*diag((bt'*Wk*bt));
%                    ExplogLD = log(ld)+0.5*(-ld.*(1-ld)).*diag(bt'*Wk*bt);
%                    sumPPll=sumPPll+sum(dN(:,k).*ExplogLD - ExplambdaDelta); 
%                     
%                 end
% 
%                 
%             end
% 
%             logll = -Dx*K/2*log(2*pi)-K/2*log(det(Q)) ...
%                     - Dx/2*log(2*pi) -1/2*log(det(Px0))  ...
%                     +sumPPll - 1/2*trace((eye(size(Q))/Q)*sumXkTerms) ...
%                     -Dx/2;
%                 string0 = ['logll: ' num2str(logll)];
%                 disp(string0);
%                 if(DEBUG==1)
%                     string1 = ['-K/2*log(det(Q)):' num2str(-K/2*log(det(Q)))];
%                     string12= ['Constants: ' num2str(-Dx*K/2*log(2*pi)-Dx/2*log(2*pi) -Dx/2 -1/2*log(det(Px0)))];
%                     string2 = ['SumPPll: ' num2str(sumPPll)];
%                     string3 = ['-.5*trace(Q\sumXkTerms): ' num2str(-.5*trace(Q\sumXkTerms))];
%                     
%                     disp(string1);
%                     disp(['Q=' num2str(diag(Q)')]);
%                     disp(string12);
%                     disp(string2);
%                     disp(string3);
%                  
%                 end
% 
%                 ExpectationSums.Sxkm1xkm1=Sxkm1xkm1;
%                 ExpectationSums.Sxkm1xk=Sxkm1xk;
%                 ExpectationSums.Sxkxkm1=Sxkxkm1;
%                 ExpectationSums.Sxkxk=Sxkxk;
%                 ExpectationSums.sumXkTerms=sumXkTerms;
%                 ExpectationSums.sumPPll=sumPPll;
%                 ExpectationSums.Sx0 = Ex0Gy;
%                 ExpectationSums.Sx0x0 = Px0Gy + Ex0Gy*Ex0Gy';
%                 ExpectationSums.A = A;
%                 ExpectationSums.Q = Q;
%                 ExpectationSums.mu = mu;
%                 ExpectationSums.beta = beta;
%                 ExpectationSums.gamma = gamma;
% 
%         end
        function [Ahat, Qhat, muhat_new, betahat_new, gammahat_new, x0hat, Px0hat] = PP_MStep(dN, x_K,W_K,x0, Px0, ExpectationSums,fitType, muhat, betahat,gammahat, windowTimes, HkAll,PPEM_Constraints,MstepMethod)
            if(nargin<14 || isempty(MstepMethod))
                MstepMethod = 'GLM'; %GLM or NewtonRaphson
            end
            if(nargin<13 || isempty(PPEM_Constraints))
                PPEM_Constraints = DecodingAlgorithms.PP_EMCreateConstraints;
            end
           
            Sxkm1xkm1=ExpectationSums.Sxkm1xkm1;
            Sxkxkm1=ExpectationSums.Sxkxkm1;
            sumXkTerms = ExpectationSums.sumXkTerms;
            [dx,K] = size(x_K);   
            numCells=size(dN,1);
            
            if(PPEM_Constraints.AhatDiag==1)
                I=eye(dx,dx);
                Ahat = (Sxkxkm1.*I)/(Sxkm1xkm1.*I);
            else
                Ahat = Sxkxkm1/Sxkm1xkm1;
            end
           
            
            if(PPEM_Constraints.QhatDiag==1)
                 if(PPEM_Constraints.QhatIsotropic==1)
                     Qhat=1/(dx*K)*trace(sumXkTerms)*eye(dx,dx);
                 else
                     I=eye(dx,dx);
                     Qhat=1/K*(sumXkTerms.*I);
                     Qhat = (Qhat + Qhat')/2;
                 end
             else
                 Qhat=1/K*sumXkTerms;
                 Qhat = (Qhat + Qhat')/2;
             end
            
             if(PPEM_Constraints.Estimatex0)
                x0hat = (inv(Px0)+Ahat'/Qhat*Ahat)\(Ahat'/Qhat*x_K(:,1)+Px0\x0);
            else
                x0hat = x0;
            end
             
            if(PPEM_Constraints.EstimatePx0==1)
                if(PPEM_Constraints.Px0Isotropic==1)
                   Px0hat=(trace(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0'))/(dx*K))*eye(dx,dx); 
                else
                    I=eye(dx,dx);
                    Px0hat =(x0hat*x0hat' - x0*x0hat' - x0hat*x0' +(x0*x0')).*I;
                    Px0hat = (Px0hat+Px0hat')/2;
                end
                
            else
                Px0hat =Px0;
            end
             
             betahat_new =betahat;
             gammahat_new = gammahat;
             muhat_new = muhat;
             
            %Compute the new CIF beta using the GLM
            if(strcmp(fitType,'poisson'))
                algorithm = 'GLM';
            else
                algorithm = 'BNLRCG';
            end
            
            % Estimate params via GLM
            if(strcmp(MstepMethod,'GLM'))
                clear c; close all;
                time=(0:length(x_K)-1)*.001;
                labels = cell(1,dx);
                labels2 = cell(1,dx+1);
                labels2{1} = 'vel';
                for i=1:dx
                    labels{i} = strcat('v',num2str(i));
                    labels2{i+1} = strcat('v',num2str(i));
                end
                vel = Covariate(time,x_K','vel','time','s','m/s',labels);
                baseline = Covariate(time,ones(length(time),1),'Baseline','time','s','',...
                    {'constant'});
                for i=1:size(dN,1)
                    spikeTimes = time(find(dN(i,:)==1));
                    nst{i} = nspikeTrain(spikeTimes);
                end
                nspikeColl = nstColl(nst);
                cc = CovColl({vel,baseline});
                trial = Trial(nspikeColl,cc);
                selfHist = windowTimes ; NeighborHist = []; sampleRate = 1000; 
                clear c;
                
                

                if(gammahat==0)
                    c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,[],NeighborHist); 
                else
                    c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,selfHist,NeighborHist); 
                end
                c{1}.setName('Baseline');
                cfgColl= ConfigColl(c);
                warning('OFF');

                results = Analysis.RunAnalysisForAllNeurons(trial,cfgColl,0,algorithm);
                temp = FitResSummary(results);
                tempCoeffs = squeeze(temp.getCoeffs);
                if(gammahat==0)
                    betahat(1:dx,:) = tempCoeffs(2:(dx+1),:);
                    muhat = tempCoeffs(1,:)';
                else
                    betahat(1:dx,:) = tempCoeffs(2:(dx+1),:);
                    muhat = tempCoeffs(1,:)';
                    histTemp = squeeze(temp.getHistCoeffs);
                    histTemp = reshape(histTemp, [length(windowTimes)-1 numCells]);
                    histTemp(isnan(histTemp))=0;
                    gammahat=histTemp;
                end
            else
                
            % Estimate via Newton-Raphson
                           % Estimate via Newton-Raphson
                 fprintf(['****M-step for beta**** \n']);
                 McExp=50;    
                 xKDrawExp = zeros(size(x_K,1),K,McExp);
                 diffTol = 1e-5;

                % Generate the Monte Carlo samples
                for k=1:K
                    WuTemp=(W_K(:,:,k));
                    [chol_m,p]=chol(WuTemp);
                    z=normrnd(0,1,size(x_K,1),McExp);
                    xKDrawExp(:,k,:)=repmat(x_K(:,k),[1 McExp])+(chol_m*z);
                end
                
                               % Stimulus Coefficients
                pool = matlabpool('size');
                if(pool==0)
                    for c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
                        fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
                            if(iter==1)
                                fprintf('%d',iter);
                            else
                                fprintf(',%d',iter);
                            end
                            if(strcmp(fitType,'poisson'))
                                HessianTerm = zeros(size(x_K,1),size(x_K,1));
                                GradTerm = zeros(size(x_K,1),1);
                                xkPerm = permute(xKDraw,[2 3 1]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaXkXkT = 1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    GradTerm = GradTerm+dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                    HessianTerm=HessianTerm-ExpLambdaXkXkT;

                                end

                            elseif(strcmp(fitType,'binomial'))
                                HessianTerm = zeros(size(x_K,1),size(x_K,1));
                                GradTerm = zeros(size(x_K,1),1);
                                xkPerm = permute(xKDraw,[1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaSquaredXk = 1/McExp*sum(repmat(ld.^2,[size(xk,1),1]).*xk,2);
                                    GradTerm = GradTerm+dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                    HessianTerm=HessianTerm+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;

                                end

                            end
                            if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                betahat_newTemp = betahat_new(:,c);
                            else
                                betahat_newTemp = (betahat_new(:,c)-HessianTerm\GradTerm);
                                if(any(isnan(betahat_newTemp)))
                                    betahat_newTemp = betahat_new(:,c);

                                end
                            end
                            mabsDiff = max(abs(betahat_newTemp - betahat_new(:,c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            betahat_new(:,c)=betahat_newTemp;
                            iter=iter+1;
                        end
                        fprintf('\n');              
                    end 
                else
                    HessianTerm = zeros(size(betahat,1),size(betahat,1),numCells);
                    GradTerm = zeros(size(betahat,1),numCells);
                    betahat_newTemp=betahat_new;
                    parfor c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
                        fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
                            if(iter==1)
                                fprintf('%d',iter);
                            else
                                fprintf(',%d',iter);
                            end
                            if(strcmp(fitType,'poisson'))
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaXkXkT = 1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    if(k==1)
                                        GradTerm(:,c) = dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                        HessianTerm(:,:,c)=-ExpLambdaXkXkT;
                                    else
                                        GradTerm(:,c) = GradTerm(:,c)+dN(c,k)*x_K(:,k) - ExpLambdaXk;
                                        HessianTerm(:,:,c)=HessianTerm(:,:,c)-ExpLambdaXkXkT;
                                    end

                                end

                            elseif(strcmp(fitType,'binomial'))
                                    xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk=xKDrawExp(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExplambdaDeltaXkXk=1/McExp*(repmat(ld,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaSqXkXkT=1/McExp*(repmat(ld.^2,[size(xk,1),1]).*xk)*xk';
                                    ExplambdaDeltaCubeXkXkT=1/McExp*(repmat(ld.^3,[size(xk,1),1]).*xk)*xk';
                                    ExpLambdaXk = 1/McExp*sum(repmat(ld,[size(xk,1),1]).*xk,2);
                                    ExpLambdaSquaredXk = 1/McExp*sum(repmat(ld.^2,[size(xk,1),1]).*xk,2);
                                    if(k==1)
                                        GradTerm(:,c) = dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                        HessianTerm(:,:,c)=ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;
                                    else
                                        GradTerm(:,c) = GradTerm(:,c)+dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExpLambdaXk+ExpLambdaSquaredXk;
                                        HessianTerm(:,:,c)=HessianTerm(:,:,c)+ExplambdaDeltaXkXk+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubeXkXkT;
                                    end
                                end

                            end
                            if(any(any(isnan(HessianTerm(:,:,c)))) || any(any(isinf(HessianTerm(:,:,c)))))
                                betahat_newTemp = betahat_new(:,c);
                            else
                                betahat_newTemp = (betahat_new(:,c)-HessianTerm(:,:,c)\GradTerm(:,c));
                                if(any(isnan(betahat_newTemp)))
                                    betahat_newTemp = betahat_new(:,c);

                                end
                            end
                            mabsDiff = max(abs(betahat_newTemp - betahat_new(:,c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            betahat_new(:,c)=betahat_newTemp;
                            iter=iter+1;
                        end
                        fprintf('\n');              
                    end 
                end
                clear GradTerm HessianTerm;
                 %Compute the CIF means 
                 if(pool==0)
                     for c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
    %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
    %                         if(iter==1)
    %                             fprintf('%d',iter);
    %                         else
    %                             fprintf(',%d',iter);
    %                         end
                            if(strcmp(fitType,'poisson'))
                                HessianTerm = zeros(size(1,1),size(1,1));
                                GradTerm = zeros(size(1,1),1);
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaDelta = 1/McExp*sum(ld,2);
                                    GradTerm = GradTerm+(dN(c,k) - ExpLambdaDelta);
                                    HessianTerm=HessianTerm-ExpLambdaDelta;

                                end

                            elseif(strcmp(fitType,'binomial'))
                                HessianTerm = zeros(size(1,1),size(1,1));
                                GradTerm = zeros(size(1,1),1);
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExpLambdaDelta =1/McExp*(sum(ld,2));
                                    ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                    ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                    GradTerm = GradTerm+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                    HessianTerm=HessianTerm+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);

                                end

                            end
                            if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                muhat_newTemp = muhat_new(c);
                            else
                                muhat_newTemp = (muhat_new(c)-HessianTerm\GradTerm);
                                if(any(isnan(muhat_newTemp)))
                                    muhat_newTemp = muhat_new(c);

                                end
                            end
                            mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            muhat_new(c)=muhat_newTemp;
                            iter=iter+1;
                        end
    %                     fprintf('\n');              
                     end 
                 else
                    HessianTerm = zeros(1,numCells);
                    GradTerm = zeros(1,numCells);
                    parfor c=1:numCells
                        converged=0;
                        iter = 1;
                        maxIter=100;
    %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                        while(~converged && iter<maxIter)
    %                         if(iter==1)
    %                             fprintf('%d',iter);
    %                         else
    %                             fprintf(',%d',iter);
    %                         end
                            if(strcmp(fitType,'poisson'))
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms);
                                    ExpLambdaDelta = 1/McExp*sum(ld,2);
                                    if(k==1)
                                        GradTerm(c) = (dN(c,k) - ExpLambdaDelta);
                                        HessianTerm(c)=-ExpLambdaDelta;
                                    else
                                        GradTerm(c) = GradTerm(c)+(dN(c,k) - ExpLambdaDelta);
                                        HessianTerm(c)=HessianTerm(c)-ExpLambdaDelta;
                                    end

                                end

                            elseif(strcmp(fitType,'binomial'))
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                for k=1:K
                                    Hk = (HkAll(:,:,c));
                                    Wk = W_K(:,:,k);
%                                     xk = squeeze(xKDrawExp(:,k,:));
                                    xk = xkPerm(:,:,k);
                                   if(size(Hk,1)==numCells)
                                       Hk = Hk';
                                   end

                                    if(numel(gammahat)==1)
                                        gammaC=gammahat;
                                    %                             gammaC=repmat(gammaC,[1 numCells]);
                                    else 
                                        gammaC=gammahat(:,c);
                                    end

                                    terms =muhat_new(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                    ld=exp(terms)./(1+exp(terms));
                                    ExpLambdaDelta =1/McExp*(sum(ld,2));
                                    ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                    ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                    if(k==1)
                                        GradTerm(c) = (dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                        HessianTerm(c)=(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);
                                    else
                                         GradTerm(c) = GradTerm(c)+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq);
                                        HessianTerm(c)=HessianTerm(c)+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed);
                                    end

                                end

                            end
                            if(any(any(isnan(HessianTerm(c)))) || any(any(isinf(HessianTerm(c)))))
                                muhat_newTemp = muhat_new(c);
                            else
                                muhat_newTemp = (muhat_new(c)-HessianTerm(c)\GradTerm(c));
                                if(any(isnan(muhat_newTemp)))
                                    muhat_newTemp = muhat_new(c);

                                end
                            end
                            mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
                            if(mabsDiff<diffTol)
                                converged=1;
                            end
                            muhat_new(c)=muhat_newTemp;
                            iter=iter+1;
                        end
    %                     fprintf('\n');              
                     end 
                 end
                 clear HessianTerm GradTerm;
                 
                 
                 %Compute the history coeffs
                 if(~isempty(windowTimes) && any(any(gammahat_new~=0)))
                     if(pool==0)
                         for c=1:numCells
                            converged=0;
                            iter = 1;
                            maxIter=100;
        %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                            while(~converged && iter<maxIter)
        %                         if(iter==1)
        %                             fprintf('%d',iter);
        %                         else
        %                             fprintf(',%d',iter);
        %                         end
                                if(strcmp(fitType,'poisson'))
                                    HessianTerm = zeros(size(gammahat,1),size(gammahat,1));
                                    GradTerm = zeros(size(gammahat,1),1);
                                    xkPerm = permute(xKDrawExp, [1 3 2]);
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end

                                        if(numel(gammahat_new)==1)
                                            gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                                        else 
                                            gammaC=gammahat_new(:,c);
                                        end

                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms);
                                        ExpLambdaDelta = 1/McExp*sum(ld,2);
                                        GradTerm = GradTerm+(dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                        HessianTerm=HessianTerm-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);

                                    end

                                elseif(strcmp(fitType,'binomial'))
                                    HessianTerm = zeros(size(gammahat,1),size(gammahat,1));
                                    GradTerm = zeros(size(gammahat,1),1);
                                    xkPerm = permute(xKDrawExp, [1 3 2]);
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk=xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end

                                        if(numel(gammahat_new)==1)
                                            gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                                        else 
                                            gammaC=gammahat_new(:,c);
                                        end

                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms)./(1+exp(terms));
                                        ExpLambdaDelta =1/McExp*(sum(ld,2));
                                        ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                        ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                        GradTerm = GradTerm+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                        HessianTerm=HessianTerm+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);

                                    end

                                end
                                if(any(any(isnan(HessianTerm))) || any(any(isinf(HessianTerm))))
                                    gammahat_newTemp = gammahat_new(:,c);
                                else
                                    gammahat_newTemp = (gammahat_new(:,c)-HessianTerm\GradTerm);
                                    if(any(isnan(gammahat_newTemp)))
                                        gammahat_newTemp = gammahat_new(:,c);

                                    end
                                end
                                mabsDiff = max(abs(gammahat_newTemp - gammahat_new(:,c)));
                                if(mabsDiff<diffTol)
                                    converged=1;
                                end
                                gammahat_new(:,c)=gammahat_newTemp;
                                iter=iter+1;
                            end
        %                     fprintf('\n');              
                         end 
                     else
                         HessianTerm = zeros(size(gammahat,1),size(gammahat,1),numCells);
                         GradTerm = zeros(size(gammahat,1),numCells);
                         parfor c=1:numCells
                            converged=0;
                            iter = 1;
                            maxIter=100;
        %                     fprintf(['neuron:' num2str(c) ' iter: ']);
                            if(numel(gammahat_new)==1)
                                gammaC=gammahat_new;
                                        %                             gammaC=repmat(gammaC,[1 numCells]);
                            else 
                                gammaC=gammahat_new(:,c);
                            end
                            while(~converged && iter<maxIter)
        %                         if(iter==1)
        %                             fprintf('%d',iter);
        %                         else
        %                             fprintf(',%d',iter);
        %                         end
                                xkPerm = permute(xKDrawExp, [1 3 2]);
                                if(strcmp(fitType,'poisson'))
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end



                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms);
                                        ExpLambdaDelta = 1/McExp*sum(ld,2);
                                        if(k==1)
                                            GradTerm(:,c) = (dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                            HessianTerm(:,:,c)=-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);
                                        else
                                            GradTerm(:,c) = GradTerm(:,c)+(dN(c,k) - ExpLambdaDelta)*Hk(k,:)';
                                            HessianTerm(:,:,c)=HessianTerm(:,:,c)-ExpLambdaDelta*Hk(k,:)'*Hk(k,:);
                                        end
                                    end

                                elseif(strcmp(fitType,'binomial'))
                                    for k=1:K
                                        Hk = (HkAll(:,:,c));
                                        Wk = W_K(:,:,k);
%                                         xk = squeeze(xKDrawExp(:,k,:));
                                        xk = xkPerm(:,:,k);
                                       if(size(Hk,1)==numCells)
                                           Hk = Hk';
                                       end


                                        terms =muhat(c)+betahat(:,c)'*xk+gammaC'*Hk(k,:)';
                                        ld=exp(terms)./(1+exp(terms));
                                        ExpLambdaDelta =1/McExp*(sum(ld,2));
                                        ExpLambdaDeltaSq = 1/McExp*(sum(ld.^2,2));
                                        ExpLambdaDeltaCubed = 1/McExp*(sum(ld.^3,2));
                                        if(k==1)
                                            GradTerm(:,c) = (dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                            HessianTerm(:,:,c)=(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);
                                        else
                                            GradTerm(:,c) = GradTerm(:,c)+(dN(c,k)-(dN(c,k)+1)*ExpLambdaDelta+ExpLambdaDeltaSq)*Hk(k,:)';
                                            HessianTerm(:,:,c)=HessianTerm(:,:,c)+(-ExpLambdaDelta*(dN(c,k)+1)+ExpLambdaDeltaSq*(dN(c,k)+3)-2*ExpLambdaDeltaCubed)*Hk(k,:)'*Hk(k,:);
                                        end

                                    end

                                end
                                if(any(any(isnan(HessianTerm(:,:,c)))) || any(any(isinf(HessianTerm(:,:,c)))))
                                    gammahat_newTemp = gammaC;
                                else
                                    gammahat_newTemp = (gammaC-HessianTerm(:,:,c)\GradTerm(:,c));
                                    if(any(isnan(gammahat_newTemp)))
                                        gammahat_newTemp = gammaC;

                                    end
                                end
                                mabsDiff = max(abs(gammahat_newTemp - gammaC));
                                if(mabsDiff<diffTol)
                                    converged=1;
                                end
                                gammaC=gammahat_newTemp;
                                iter=iter+1;
                            end
                            gamma_new(:,c) =gammaC;
        %                     fprintf('\n');              
                         end 
                     end
                 end
                 clear HessianTerm GradTerm; 
            end
        end

%         function [Ahat, Qhat, muhat_new, betahat_new, gammahat_new, x0hat, Px0hat] = PP_MStep(dN,x_K,W_K,x0, ExpectationSums,fitType, muhat, betahat,gammahat, windowTimes, HkAll,MstepMethod)
%             if(nargin<12 || isempty(MstepMethod))
%                 MstepMethod = 'GLM'; %GLM or NewtonRaphson
%             end
%             Sxkm1xkm1=ExpectationSums.Sxkm1xkm1;
%             Sxkxkm1=ExpectationSums.Sxkxkm1;
%             Sxkxk=ExpectationSums.Sxkxk;
%             sumXkTerms = ExpectationSums.sumXkTerms;
%             Sx0 = ExpectationSums.Sx0;
%             Sx0x0 = ExpectationSums.Sx0x0;
%             K = size(x_K,2);   
%             numCells=size(dN,1);
%             numStates = size(x_K,1);
%             Ahat = Sxkxkm1/Sxkm1xkm1;
%         
%             Px0hat =(Sx0x0 - x0*Sx0' - Sx0*x0' +(x0*x0'));
%              
% %              [V,D] = eig(Px0hat);
% %              D(D<0)=1e-9;
% %              Px0hat = V*D*V';
%             Px0hat = (Px0hat+Px0hat')/2;
% %              Px0hat = diag(diag(Px0hat));
%             x0hat  = Sx0;
%         
%             Qhat=1/K*sumXkTerms;
% %             [V,D] = eig(Qhat);
% %             D(D<=0)=1e-9;
% %             Qhat = V*D*V';
%             Qhat = (Qhat + Qhat')/2;
%             if(det(Qhat)<=0)
%                 Qhat = ExpectationSums.Q; % Keep prior value
%             end
%             
%                
%              
%             betahat_new =betahat;
%             gammahat_new = gammahat;
%             muhat_new = muhat;
%              
%             %Compute the new CIF beta using the GLM
%             if(strcmp(fitType,'poisson'))
%                 algorithm = 'GLM';
%             else
%                 algorithm = 'BNLRCG';
%             end
%             
%             % Estimate params via GLM
%             if(strcmp(MstepMethod,'GLM'))
%                 clear c; close all;
%                 time=(0:length(x_K)-1)*.001;
%                 labels = cell(1,numStates);
%                 labels2 = cell(1,numStates+1);
%                 labels2{1} = 'vel';
%                 for i=1:numStates
%                     labels{i} = strcat('v',num2str(i));
%                     labels2{i+1} = strcat('v',num2str(i));
%                 end
%                 vel = Covariate(time,x_K','vel','time','s','m/s',labels);
%                 baseline = Covariate(time,ones(length(time),1),'Baseline','time','s','',...
%                     {'constant'});
%                 for i=1:size(dN,1)
%                     spikeTimes = time(dN(i,:)==1);
%                     nst{i} = nspikeTrain(spikeTimes);
%                 end
%                 nspikeColl = nstColl(nst);
%                 cc = CovColl({vel,baseline});
%                 trial = Trial(nspikeColl,cc);
%                 selfHist = windowTimes ; NeighborHist = []; sampleRate = 1000; 
%                 clear c;
%                 
%                 
% 
%                 if(gammahat==0)
%                     c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,[],NeighborHist); 
%                 else
%                     c{1} = TrialConfig({{'Baseline','constant'},labels2},sampleRate,selfHist,NeighborHist); 
%                 end
%                 c{1}.setName('Baseline');
%                 cfgColl= ConfigColl(c);
%                 warning('OFF');
% 
%                 results = Analysis.RunAnalysisForAllNeurons(trial,cfgColl,0,algorithm);
%                 temp = FitResSummary(results);
%                 tempCoeffs = squeeze(temp.getCoeffs);
%                 if(gammahat==0)
%                     betahat(1:numStates,:) = tempCoeffs(2:(numStates+1),:);
%                     muhat = tempCoeffs(1,:)';
%                 else
%                     betahat(1:numStates,:) = tempCoeffs(2:(numStates+1),:);
%                     muhat = tempCoeffs(1,:)';
%                     histTemp = squeeze(temp.getHistCoeffs);
%                     histTemp = reshape(histTemp, [length(windowTimes)-1 numCells]);
%                     histTemp(isnan(histTemp))=0;
%                     gammahat=histTemp;
%                     if(size(gammahat,2)~=size(muhat,1))
%                         gammahat = gammahat';
%                     end
%                 end
%             else
%                 
%             % Estimate via Newton-Raphson
%                  fprintf(['****M-step for beta**** \n']);
%                  for c=1:numCells
%     %                  c
% 
% 
%                      converged=0;
%                      iter = 1;
%                      maxIter=100;
%     %                  disp(['M-step for beta, neuron:' num2str(c) ' iter: ' num2str(c) ' of ' num2str(maxIter)]); 
%                      fprintf(['neuron:' num2str(c) ' iter: ']);
%                      while(~converged && iter<maxIter)
% 
%                         if(iter==1)
%                             fprintf('%d',iter);
%                         else
%                             fprintf(',%d',iter);
%                         end
%                         if(strcmp(fitType,'poisson'))
%                             gradQ=zeros(size(betahat_new(:,c),1),1);
%                             jacQ =zeros(size(betahat_new(:,c),1),size(betahat_new(:,c),1));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 xk = x_K(:,k);
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
%                                 ld=exp(terms);
% 
%                                 numStates =length(xk);
%                                 ExplambdaDeltaXk = zeros(numStates,1);
%                                 ExplambdaDeltaXkXkT = zeros(numStates,numStates);
%                                 for m=1:numStates
%                                      sm = zeros(numStates,1);
%                                      sm(m) =1;
%                                      bt=betahat_new(:,c);
%                                      ExplambdaDeltaXk(m) = ld*sm'*xk+...
%                                          .5*trace(ld*(bt*xk'*sm*bt'+sm*bt'+bt*sm')*Wk);
%                                     for n=1:m
%                                         sn = zeros(numStates,1);
%                                         sn(n) =1; 
%                                         ExplambdaDeltaXkXkT(n,m) = ld*xk'*sm*sn'*xk+...
%                                             +trace(ld*(2*bt*xk'*sn*sm'*xk*bt'+bt*xk'*sn*sm'+sn*sm'*xk*bt'+sn*sm')*Wk);
%                                         if(n~=m)
%                                             ExplambdaDeltaXkXkT(n,m)=ExplambdaDeltaXkXkT(m,n);
%                                         end
%                                     end
%                                 end
% 
%                                 gradQ = gradQ + (dN(c,k)*xk - ExplambdaDeltaXk);
%                                 jacQ  = jacQ  - ExplambdaDeltaXkXkT;
%                             end
% 
% 
%                         elseif(strcmp(fitType,'binomial'))
%                             gradQ=zeros(size(betahat_new(:,c),1),1);
%                             jacQ =zeros(size(betahat_new(:,c),1),size(betahat_new(:,c),1));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 xk = x_K(:,k);                    
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms =muhat(c)+betahat_new(:,c)'*xk+gammaC'*Hk(k,:)';
%                                 ld=exp(terms)./(1+exp(terms));
% 
%                                 numStates =length(xk);
%                                 ExplambdaDeltaXk = zeros(numStates,1);
%                                 ExplambdaDeltaSqXk = zeros(numStates,1);
%                                 ExplambdaDeltaXkXkT = zeros(numStates,numStates);
%                                 ExplambdaDeltaSqXkXkT = zeros(numStates,numStates);
%                                 ExplambdaDeltaCubedXkXkT = zeros(numStates,numStates);
%                                 for m=1:numStates
%                                      sm = zeros(numStates,1);
%                                      sm(m) =1;
%                                      bt=betahat_new(:,c);
%                                      ExplambdaDeltaXk(m) = ld*sm'*xk+...
%                                          +.5*trace(ld*(bt*xk'*sm*bt'+sm*bt'+bt*sm')*Wk)...
%                                          -.5*trace((ld^2)*(3*bt*xk'*sm*bt'+sm*bt'+bt*sm')*Wk)...
%                                          +.5*trace((ld^3)*(2*bt*xk'*sm*bt')*Wk);
%                                      ExplambdaDeltaSqXk(m) = (ld)^2*sm'*xk+...
%                                          +trace((ld^2)*(2*bt*xk'*sm*bt'+sm*bt'+bt*sm')*Wk)...
%                                          -trace((ld^3)*(2*bt*xk'*sm*bt'+3*bt*xk'*sm*bt'+sm*bt'+bt*sm')*Wk)...
%                                          +trace(3*(ld^4)*(bt*xk'*sm*bt')*Wk);
% 
%                                     for n=1:m
%                                         sn = zeros(numStates,1);
%                                         sn(n) =1; 
%                                         ExplambdaDeltaXkXkT(n,m) = ld*xk'*sm*sn'*xk+...
%                                             +0.5*trace((ld)*(bt*xk'*sn*sm'*xk*bt'+2*sn*sm'*xk*bt'+2*bt*xk'*sn*sm'+2*sn*sm')*Wk)...
%                                             -0.5*trace((ld)^2*(3*bt*xk'*sn*sm'*xk*bt'+2*sn*sm'*xk*bt'+2*bt*xk'*sn*sm')*Wk)...
%                                             +0.5*trace((ld)^3*(2*bt*xk'*sn*sm'*xk*bt')*Wk);
%                                         ExplambdaDeltaSqXkXkT(n,m) = (ld)^2*xk'*sm*sn'*xk+...
%                                             +trace((ld)^2*(2*bt*xk'*sn*sm'*xk*bt'+2*sn*sm'*xk*bt'+2*bt*xk'*sn*sm'+sn*sm')*Wk)...
%                                             -trace((ld)^3*(5*bt*xk'*sn*sm'*xk*bt'+2*sn*sm'*xk*bt'+2*bt*xk'*sn*sm')*Wk)...
%                                             +trace((ld)^4*(3*bt*xk'*sn*sm'*xk*bt')*Wk);
% 
%                                         ExplambdaDeltaCubedXkXkT(n,m) = (ld)^3*xk'*sm*sn'*xk+...
%                                             +0.5*trace((ld)^3*(9*bt*xk'*sn*sm'*xk*bt'+6*sn*sm'*xk*bt'+6*bt*xk'*sn*sm'+2*sn*sm')*Wk)...
%                                             -0.5*trace((ld)^4*(21*bt*xk'*sn*sm'*xk*bt'+6*sn*sm'*xk*bt'+6*bt*xk'*sn*sm')*Wk)...
%                                             +0.5*trace((ld)^5*(12*bt*xk'*sn*sm'*xk*bt')*Wk);
% 
%                                         if(n~=m)
%                                             ExplambdaDeltaXkXkT(n,m)=ExplambdaDeltaXkXkT(m,n);
%                                             ExplambdaDeltaSqXkXkT(n,m)=ExplambdaDeltaSqXkXkT(m,n);
%                                             ExplambdaDeltaCubedXkXkT(n,m)=ExplambdaDeltaCubedXkXkT(m,n);
%                                         end
%                                     end
%                                 end
% 
%                                 gradQ = gradQ + dN(c,k)*x_K(:,k) - (dN(c,k)+1)*ExplambdaDeltaXk+ExplambdaDeltaSqXk;
%                                 jacQ  = jacQ  + ExplambdaDeltaXkXkT+ExplambdaDeltaSqXkXkT-2*ExplambdaDeltaCubedXkXkT;
%                             end
%                         end
% 
% 
%     %                    gradQ=0.01*gradQ;
% 
% 
%                         if(any(any(isnan(jacQ))) || any(any(isinf(jacQ))))
%                             betahat_newTemp = betahat_new(:,c);
%                         else
%                             betahat_newTemp = (betahat_new(:,c)-jacQ\gradQ);
%                             if(any(isnan(betahat_newTemp)))
%                                 betahat_newTemp = betahat_new(:,c);
% 
%                             end
%                         end
%                         mabsDiff = max(abs(betahat_newTemp - betahat_new(:,c)));
%                         if(mabsDiff<10^-2)
%                             converged=1;
%                         end
%                         betahat_new(:,c)=betahat_newTemp;
%                         iter=iter+1;
%                     end
%                     fprintf('\n');              
%                  end 
% 
% 
%                  %Compute the new CIF means
%                  muhat_new =muhat;
%                  for c=1:numCells
%                      converged=0;
%                      iter = 1;
%                      maxIter=100;
%                      while(~converged && iter<maxIter)
%                         if(strcmp(fitType,'poisson'))
%                             gradQ=zeros(size(muhat_new(c),2),1);
%                             jacQ =zeros(size(muhat_new(c),2),size(muhat_new(c),2));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                 ld = exp(terms);
%                                 bt = betahat(:,c);
%                                 ExplambdaDelta =ld +0.5*trace(ld*bt*bt'*Wk);
% 
% 
%                                 gradQ = gradQ + dN(c,k)' - ExplambdaDelta;
%                                 jacQ  = jacQ  - ExplambdaDelta;
%                             end
% 
% 
%                         elseif(strcmp(fitType,'binomial'))
%                             gradQ=zeros(size(muhat_new(c),2),1);
%                             jacQ =zeros(size(muhat_new(c),2),size(muhat_new(c),2));
%                             for k=1:K
% %                                 Hk=HkAll{c};
%                                 Hk = squeeze(HkAll(:,:,c));
%                                 Wk = W_K(:,:,k);
%                                 if(numel(gammahat)==1)
%                                     gammaC=gammahat;
%                                 else 
%                                     gammaC=gammahat(:,c);
%                                 end
%                                 terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                 ld = exp(terms)./(1+exp(terms));
%                                 bt = betahat(:,c);
%                                 ExplambdaDelta = ld+0.5*trace(bt*bt'*(ld)*(1-ld)*(1-2*ld)*Wk);
%                                 ExplambdaDeltaSq = (ld)^2+...
%                                     0.5*trace((ld)^2*(1-ld)*(2-3*ld)*bt*bt'*Wk);
%                                 ExplambdaDeltaCubed = (ld)^3+...
%                                     0.5*trace(3*(ld)^3*(3-7*ld+4*(ld)^2)*bt*bt'*Wk);
% 
%                                 gradQ = gradQ + dN(c,k)' -(dN(c,k)+1)*ExplambdaDelta...
%                                     +ExplambdaDeltaSq;
%                                 jacQ  = jacQ  - (dN(c,k)+1)*ExplambdaDelta...
%                                     +(dN(c,k)+3)*ExplambdaDeltaSq...
%                                     -3*ExplambdaDeltaCubed;
%                             end
% 
%                         end
%     %                     gradQ=0.01*gradQ;
%                         muhat_newTemp = (muhat_new(c)'-(1/jacQ)*gradQ)';
%                         if(any(isnan(muhat_newTemp)))
%                             muhat_newTemp = muhat_new(c);
% 
%                         end
%                         mabsDiff = max(abs(muhat_newTemp - muhat_new(c)));
%                         if(mabsDiff<10^-2)
%                             converged=1;
%                         end
%                         muhat_new(c)=muhat_newTemp;
%                         iter=iter+1;
%                      end
% 
%                 end
% 
%     %             Compute the history parameters
%                 gammahat_new = gammahat;
%                 if(~isempty(windowTimes) && any(any(gammahat_new~=0)))
%                      for c=1:numCells
%                          converged=0;
%                          iter = 1;
%                          maxIter=100;
%                          while(~converged && iter<maxIter)
%                             if(strcmp(fitType,'poisson'))
%                                 gradQ=zeros(size(gammahat_new(c),2),1);
%                                 jacQ =zeros(size(gammahat_new(c),2),size(gammahat_new(c),2));
%                                 for k=1:K
% %                                     Hk=HkAll{c};
%                                     Hk = squeeze(HkAll(:,:,c));
%                                     Wk = W_K(:,:,k);
%                                     if(numel(gammahat)==1)
%                                         gammaC=gammahat;
%                                     else 
%                                         gammaC=gammahat(:,c);
%                                     end
%                                     terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                     ld = exp(terms);
%                                     bt = betahat(:,c);
%                                     ExplambdaDelta =ld +0.5*trace(bt*bt'*ld*Wk);
% 
% 
%                                     gradQ = gradQ + (dN(c,k)' - ExplambdaDelta)*Hk;
%                                     jacQ  = jacQ  - ExplambdaDelta*Hk*Hk';
%                                 end
% 
% 
%                             elseif(strcmp(fitType,'binomial'))
%                                 gradQ=zeros(size(gammahat_new(c),2),1);
%                                 jacQ =zeros(size(gammahat_new(c),2),size(gammahat_new(c),2));
%                                 for k=1:K
% %                                     Hk=HkAll{c};
%                                     Hk = squeeze(HkAll(:,:,c));
%                                     Wk = W_K(:,:,k);
%                                     if(numel(gammahat)==1)
%                                         gammaC=gammahat;
%                                     else 
%                                         gammaC=gammahat(:,c);
%                                     end
%                                     terms=muhat_new(c)+betahat(:,c)'*x_K(:,k)+gammaC'*Hk(k,:)';
%                                     ld = exp(terms)./(1+exp(terms));
%                                     bt = betahat(:,c);
%                                     ExplambdaDelta =ld...
%                                         +0.5*trace(bt*bt'*ld*(1-ld)*(1-2*ld)*Wk);
%                                     ExplambdaDeltaSq=ld^2 ...
%                                         +trace((ld^2*(1-ld)*(2-3*ld)*bt*bt')*Wk);
%                                     ExplambdaDeltaCubed=ld^3 ...
%                                         +0.5*trace((9*(ld^3)*(1-ld)^2*bt*bt'-3*(ld^4)*(1-ld)*bt*bt')*Wk);
%                                     gradQ = gradQ + (dN(c,k) - (dN(c,k)+1)*ExplambdaDelta+ExplambdaDeltaSq)*Hk;
%                                     jacQ  = jacQ  + -ExplambdaDelta*(dN(c,k)+1)*Hk*Hk'...
%                                         +ExplambdaDeltaSq*(dN(c,k)+3)*Hk*Hk'...
%                                         -ExplambdaDeltaCubed*2*Hk*Hk';
%                                 end
% 
%                             end
% 
% 
%     %                         gradQ=0.01*gradQ;
% 
%                             gammahat_newTemp = (gammahat_new(:,c)-(eye(size(Hk,2),size(Hk,2))/jacQ)*gradQ');
%                             if(any(isnan(gammahat_newTemp)))
%                                 gammahat_newTemp = gammahat_new(:,c);
% 
%                             end
%                             mabsDiff = max(abs(gammahat_newTemp - gammahat_new(:,c)));
%                             if(mabsDiff<10^-2)
%                                 converged=1;
%                             end
%                             gammahat_new(:,c)=gammahat_newTemp;
%                             iter=iter+1;
%                          end
% 
%                     end
%     %                  gammahat(:,c) = gammahat_new;
%                 end
% %              betahat =betahat_new;
% %              gammahat = gammahat_new;
% %              muhat = muhat_new;
%             end           
%         end
    end
end

