
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMATRS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = IMATRS(fun,x,low,upp,st,tune);
%
% solve unconstrained noisy black box optimization problem 
%    min f(x) 
%  
% fun      % function handle for f(.)
% x        % starting point (must be specified)
% low      % lower bound 
% upp      % upper bound
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax (default: inf)
%  .nfmax        %   stop if nf>=nfmax   (default: inf)
%  .qf           %   stop if qf<=accf    (default: 1e-4)
%  .fbest        %   optimal function value (default: 0)
%  .prt          %   printlevel (default: 0)
%                    0: nothing, 1: litte, >=2: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% xbest    % best point found 
% fbest    % function value at best point found 

function [xbest,fbest] = IMATRS(fun,x,tune)
global feinfo
if nargin<2
  message = 'fun,x, should be as input';
  feinfo.error= message; 
  return
end
% check function handle
if isempty(fun)
  message = 'subTR needs the function handle fun to be defined';
  feinfo.error= message; 
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  feinfo.error= message;
  return
end
% starting point
if isempty(x)
  message = 'starting point must be defined';
  feinfo.error= message;
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  feinfo.error= message;
  return       
end
% dimension
n = length(x);


if length(feinfo.upp)~=n
  nlow=n,nupp=length(feinfo.upp)
  feinfo.error='dimension mismatch'; 
end
if n==0
  % no variables
  xbest=zeros(0);fbest=inf;
  feinfo.error='no variables';
  return;
end
if ~all(feinfo.low<=feinfo.upp) 
  % feasible domain empty
  xbest=NaN;fbest=NaN;
  feinfo.error='feasible domain empty';
  return;
end
% artificial lower and upper bounds to help prevent overflow
% initialize structure containing all tuning parameters 
%
% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input

if ~exist('tune'), tune=[]; end

if ~isfield(tune,'lambda'), tune.lambda = n; end
if ~isfield(tune,'mu'), tune.mu = floor(tune.lambda/2); end
if ~isfield(tune,'sigma'), tune.sigma = 1; end

if ~isfield(tune,'zeta'), tune.zeta = 1e-20; end
if ~isfield(tune,'nu'), tune.nu = 2; end
if ~isfield(tune,'theta'), tune.theta = 2; end
if ~isfield(tune,'Deltainit'), tune.Deltainit = 10;  end
if ~isfield(tune,'Deltabar'), tune.Deltabar = 3;  end
if ~isfield(tune,'sigmamax'), tune.sigmamax = 1e10;  end

% factor for adjusting Y
if ~isfield(tune,'gammaX'), tune.gammaX = 1e3; end
% factor for adjusting gradient
if ~isfield(tune,'gammav'), tune.gammav = 1e2; end

lambda = tune.lambda; sigma=tune.sigma;

echi     = sqrt(n)*(1 - 1/n/4 - 1/n/n/21);
Parent.y = x;
s = zeros(n, 1);
M = eye(n);
ordering = 'ascend';
xf=Inf*ones(feinfo.nfmax,n+1);
eta         = 1.5;
Phalton     = haltonset(n);
ihalton     = 7;
allones     = 0;
x = max(feinfo.low,min(feinfo.upp,x));
xbest = x;
[f] = fun(x); 
fbest=f;
if feinfo.done; return; end 
nx=size(x,1);
xf(feinfo.nf,1:nx)=x';
xf(feinfo.nf,nx+1)=fbest;
% check stopping test
a = ones(lambda,1); succ = zeros(1,n); D=eye(n,n); 
it = 0; aold=inf; prt = feinfo.prt;
while 1
       it = it+1;
       x=xbest; 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % phase I: Perform integer mutation
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if prt>=2, fprintf(['start of ',num2str(it),'th inMutation\n'])
       end
       [xbest,fbest,xf,a,good,allones,succ,OffspringPop,restart]= ...
           intMutation(fun,n,xbest,fbest,a,succ,allones,D,nx,M,xf,tune);
       if prt>=2, fprintf(['end of ',num2str(it),'th inMutation\n'])
       end
       if feinfo.done; return; end 
       
        % enrich D if needed
        if (norm(xbest-x) <= 0) && (max(a) == 1) && (aold == 1)
            allones=allones+1; iexit = 0;
            while iexit == 0
                % enrich set D
                 [Phalton,ihalton,D, succ, a, iexit] = ...
                        generate_dirs(n,D,succ,a,eta,Phalton,ihalton,0);
               if iexit == 0, eta = eta + 0.5; end
               if eta >= 0.5*(norm(feinfo.upp-feinfo.low)./2),
                   iexit = 1;
               end
            end
        end
        aold = max(a);
        mu = length(OffspringPop); lambda = 2*mu;
        if restart, s = zeros(n, 1); M = eye(n); end
        esok = 0; 
       if mu>=2
         wi_raw   = log(lambda/2 + 0.5) - log((1:mu));
         wi       = wi_raw/sum(wi_raw);
         mu_eff   = max(1+eps,1/sum(wi .^2));
         c_s      = min(1.999,(mu_eff + 2)/(n + mu_eff + 5));
         c_1      = 2/((n+1.3)^2 + mu_eff);
         c_mu     = .... 
            min( 1-c_1, 2*(mu_eff - 2 + 1/mu_eff)/((n + 2)^2 + mu_eff));
         d_s      = 1 + c_s + 2*max( 0, sqrt((mu_eff-1)/(n+1)) - 1 );
         sqrt_s   = sqrt(c_s*(2-c_s)*mu_eff);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % phase II: perform selection
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if prt>=2, fprintf(['start of ',num2str(it),'th selection\n'])
          end
         ranks = RankPop(OffspringPop, ordering);
         sum_z = zeros(n, 1);
         sum_d = zeros(n, 1);
         for m = 1:mu; 
             sum_z = sum_z + wi(m)*OffspringPop(ranks(m)).sd;
             sum_d = sum_d + wi(m)*OffspringPop(ranks(m)).p;
          end
          Parent.Ps=sum_z;
           if prt>=2, fprintf(['end of ',num2str(it),'th selection\n'])
          end
          % update M
          s = (1-c_s)*s + sqrt_s*Parent.Ps;
          M = (1 - 0.5*c_1 - 0.5*c_mu) * M + (0.5*c_1)*(M*s)*s';
          for m = 1:mu
              M = M + ((0.5*c_mu*wi(m))*OffspringPop(ranks(m)).p)*...
                 OffspringPop(ranks(m)).sd';
          end  
          M=adjustY(M,tune);
         % update sigma
         pow = norm(s)/echi - 1; sigma = sigma*exp((c_s/d_s)*pow);
         sigma = round(min(tune.sigmamax,max(1,sigma)));
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % perform integer recombination
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if prt>=2
             fprintf(['start of ',num2str(it),'th intRecombination\n'])
          end
         [xbest,fbest,nx,xf,fes,yes,esok,good]=...
            intRecombination(fun,n,good,mu,sum_d,xbest,fbest,nx,xf,tune);
         if prt>=2
             fprintf(['end of ',num2str(it),'th intRecombination\n'])
         end
       end
       okTR=(((mu<=1&&~good)||(mu>1&&~good&&norm(d)==0))&&feinfo.nf>=n+1); 
       if okTR % perform integer trust region algorithm
           if prt>=2
             fprintf(['start of ',num2str(it),'th intTRS\n'])
          end
          [xbest,fbest,good,xf]=intTRS(fun,n,xbest,fbest,nx,xf,M,tune);
          if prt>=2
             fprintf(['end of ',num2str(it),'th intTRS\n'])
          end
          if feinfo.done; return; end 
      end
      if esok && ~good, xbest  = yes; fbest = fes; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
