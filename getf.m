%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f]=getf(z);
% computes function value of func at z
% collects statistics and enforces stopping test
%
% z       % argument
%
% f      % function value
%
function [f]=getf(z)

global feinfo 

% feinfo        % structure containing function info and statistics
%   .func          %   function handle
%   .low           %   lower bounds with integer variables
%   .upp           %   upper bounds with integer variables
%   .lowc          %   lower bounds with continuous variables
%   .uppc          %   upper bounds with continuous variables
%   .accf          %   desired accuracy for qf
%   .nfmax         %   bounds on number of function values
%   .secmax        %   bounds on number of seconds
%   .nf            %   number of function evaluations
%   .f             %   function value
%   .ffbest        %   best function value known to us used in qf
%   .fbest         %   best function value
%   .zfbest        %   point with best function value
%   .nfbest        %   nf of best function value
%   .zqfbest       %   point with best qf
%   .fqfbest       %   function value at best qf
%   .nfqfbest      %   nf of best qf
%   .maxInfeas     %   maximal infeasibility of evaluation points
%   .error         %   error message


feinfo.error=[]; feinfo.done=0;
if feinfo.qfbest<=feinfo.accf 
  feinfo.error='accuracy reached';
   feinfo.done = 1;
elseif feinfo.fbest<=-1e+40,
    feinfo.error='function -1e+40 reached';
    feinfo.done = 1;
elseif feinfo.nf>=feinfo.nfmax
      feinfo.error='nfmax reached';
      feinfo.done = 1;
elseif cputime-feinfo.cpu0>=feinfo.secmax, 
    feinfo.error='secmax reached'; 
    feinfo.done = 1;
end

zz = max(feinfo.low,min(z(:),feinfo.upp));
infeas=norm(zz-z,inf);
if infeas>0, [zz; z']; end
feinfo.maxInfeas=max(feinfo.maxInfeas,infeas);
y    = feinfo.clow+z.*(feinfo.cupp-feinfo.clow)./100;  
[f] = feinfo.func(y);
if isnan(f), f=inf; end;
f=max(-1e50,min(f,1e50));
fs = f; %  true function value

% add noise to the true function value
if feinfo.noisefun
   switch feinfo.noise.distr
       case 1, noise.epsilon = (2*rand-1)*feinfo.noise.level;
       case 2, noise.epsilon = randn*feinfo.noise.level;
       otherwise
        disp('error: noisy.distr should be either uniform or Gauss')
        disp('Go to CUTESTall and choose truly noisy.distr')
        next
   end
   switch feinfo.noise.type
       case 1
           f = f+noise.epsilon;
       case  2
           f= f*(1+noise.epsilon);
       otherwise 
           disp('error: noisy.type should be either ab or rel ')
           disp('Go to CUTESTall and choose truly noisy.type')
   end
end
feinfo.nf=feinfo.nf+1;

if feinfo.nf==1, feinfo.finit = fs; end;

feinfo.f=f;
feinfo.z=z;
if infeas==0 && f<=feinfo.fbest
  feinfo.fbest  = f;  
  feinfo.zbest  = z;  
  feinfo.nfbest = feinfo.nf;
  feinfo.fec=cputime-feinfo.cpu0;
  if feinfo.prt>0
    format long
    disp(['inexact function value improved at nf=',num2str(feinfo.nf),...
        ' to f=',num2str(feinfo.fbest)]) 
    format short 
  end
  finit            = feinfo.finit;
  ffbest           = max(-abs(finit)-1e12,feinfo.ffbest);
  qf               = (fs-ffbest)/((finit-ffbest)+realmin);
  qf               = max(qf,0);
  feinfo.qf        = qf;
  feinfo.zqfbest   = z;
  feinfo.qfbest    = qf;
  feinfo.fqfbest   = f;
  feinfo.nfqfbest  = feinfo.nf;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%