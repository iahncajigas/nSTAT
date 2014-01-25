function [X,yframe] = GenCovMat(t,J,y,K)

% Data parameters
trialen = 3000; %trialen = 5;
numtrials = length(t)/trialen;
avsples = max(J,K+1);
numsples = length(y);
framelength = trialen - avsples + 1;
numframes = floor((numsples - avsples + 1)/framelength);

% if J <= 0
%     disp('Output must at least depend on covariates!');
%     return;
% end
% 
% if framelength < J + K + 1
%     fprintf('\nIll-posed problem:\n');
%     fprintf('\tYou need at least as many observations as the number of parameters you are trying to estimate!\n');
%     return;
% end
% 
% 
% fprintf('\nEstimating parameters:\n');
% fprintf('\tThe length of the covariate kernel is %d and that of the output kernel %d\n',J,K);
% fprintf('\tThe data consists of a total of %d samples.\n',numsples);
% fprintf('\tThe 1st %d samples were assumed to be available.\n',avsples);
% fprintf('\tThere are a total of %d frame(s) to be processed, each of length %d.\n',numframes,framelength);


startframe = avsples;% + (i-1)*framelength;
endframe = avsples + framelength - 1;%i*framelength -1;

% Build data vector

T = reshape(t,trialen,numtrials)';
Y = reshape(y,trialen,numtrials)';

X = zeros(numtrials*framelength,J+K+1);
X(:,1) = ones(numtrials*framelength,1);
yframe = zeros(numtrials*framelength,1);

for i=1:numtrials

    for j=0:J-1
        X((i-1)*framelength + 1:i*framelength,j+K+1+1) = T(i,avsples-j:avsples-j+framelength-1)';
    end

    for k=1:K
        X((i-1)*framelength + 1:i*framelength,k+1) = Y(i,avsples-k:avsples-k+framelength-1)';
    end
    
    yframe((i-1)*framelength+1:i*framelength) = Y(i,avsples:end);
 
end

%Xframe = X(:,2:end);
X = X(:,2:end);
