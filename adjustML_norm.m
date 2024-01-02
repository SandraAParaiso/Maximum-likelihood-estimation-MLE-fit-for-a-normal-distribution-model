% Function that performs a maximum likelihood estimation (MLE) fit for a normal distribution model.
% It is used to adjust parameters (mu and sigma) of a normal distribution to observed data.
% The fit is based on comparing observed responses with predictions from the normal model.

% The code calculates the probability of responding "yes" (pYes) based on the duration of an event, using parameters mu1 and sigma1.
% Subsequently, adjustments are made to prevent issues with log(0).
% Finally, it computes the negative log-likelihood function, which is maximized during parameter fitting.

function f=adjustML_normal(P)
global simulation 
mu1=P(1);
sigma1=P(2);
pg1=0.5; 
pl1=0;
%pg=0.5; % 0.5 for 2AFC 0.0 for 1AFC -> yes-no

duration=simulation(:,1);
percentage1=simulation(:,2);
nYes=simulation(:,3);
nNo=simulation(:,4);
nTOTAL=simulation(:,5);

if pl1<0 | pl1>0.06 
   f=1e25;
   return;
end

pYes=pg1+((1-pg1-pl1).*normcdf(duration,mu1,sigma1));

 
% To eliminate the log(0) problem;
  tol = 1e-4;
  z_index = find(pYes == 0);
  if (~isempty(z_index))
    pYes(z_index) = tol*ones(length(z_index),1);
  end
  o_index = find(pYes == 1);
  if (~isempty(o_index))
    pYes(o_index) = (1-tol)*ones(length(o_index),1);
 end
 
tmp = nYes.*log(pYes) + nNo.*log(1 - pYes);
f = -sum(tmp);