%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function initf(func,low,upp,clow,cupp,st); 
% initializes global data structure for feinfo.m
% 
% func          % function handle
%   .low        %   lower bounds with integer variables
%   .upp        %   upper bounds with integer variables
%   .clow       %   lower bounds with continuous variables
%   .cupp       %   upper bounds with continuous variables
% st            % run time limits
%   .accf       %   desired accuracy for qf
%   .nfmax      %   bounds on number of function values
%   .secmax     %   bounds on number of seconds
%   .prt        %   printlevel (default: 0)
%                   0: nothing, 1: litte, >=2: more and more
%
function initf(func,low,upp,clow,cupp,st)

global feinfo
% feinfo     % structure containing function info and statistics
%   .func          %   function handle
%   .low           %   lower bounds with integer variables
%   .upp           %   upper bounds with integer variables
%   .clow          %   lower bounds with continuous variables
%   .cupp          %   upper bounds with continuous variables
%   .accf          %   desired accuracy for qf
%   .nfmax         %   bounds on number of function values
%   .secmax        %   bounds on number of seconds
%   .nf            %   number of function evaluations
%   .f             %   function value
%   .ffbest        %   best function value known to us used in qf
%   .fbest         %   best function value
%   .zfbest        %   point with best function value
%   .zqfbest       %   point with best qf
%   .fqfbest       %   function value at best qf
%   .nfqfbest      %   nf of best qf
%   .maxInfeas     %   maximal infeasibility of evaluation points
%   .error         %   error message
%   .cpu0          %   cputime at start of process

feinfo=[];
feinfo.noisefun = st.noisefun;
if feinfo.noisefun
   feinfo.noise  = st.noise;
end
feinfo.func=func;
feinfo.low=low; 
feinfo.upp=upp;
feinfo.clow=clow; 
feinfo.cupp=cupp;
feinfo.accf=st.accf;
feinfo.nfmax=st.nfmax;
feinfo.secmax=st.secmax;
feinfo.acc=NaN;
feinfo.nf=0;
feinfo.f=NaN;
feinfo.fbest=inf;
feinfo.z=nan+low;
feinfo.zbest=nan+low;
feinfo.f=NaN;
feinfo.qf=inf;
feinfo.zqfbest=NaN+low;
feinfo.fqfbest=NaN;
feinfo.nfqfbest=0;
feinfo.qfbest = NaN;
feinfo.maxInfeas=0;
feinfo.error='';
feinfo.zinit   = st.z;
feinfo.ffbest  = st.fbest;
feinfo.cpu0    = cputime;
feinfo.prt     = st.prt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
