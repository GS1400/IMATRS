%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% driverIMATRS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end

dline = ['===============================================',...
            '=============\n'];

fprintf(['IMATRS solves integer bound constrained noisy black box',...
  ' optimization problems of a not necessarily\n',...
  ' smooth function of many continuous arguments. No derivatives',...
  ' are needed. A limited amount of noise\n',...
  ' is tolerated.\n\n']);

fprintf(dline)


solverPath = input(['Insert the IMATRS path \n',...
    '>> solverPath='],'s');


if ~exist(solverPath, 'dir')
    disp('the directory does not exist')
    return
end


fprintf(dline)

global feinfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create noise

fprintf(['noise info:\n',...
         'noise.level: 0.0001/0.01/0.1 \n',...
         'noise.type:  1 (absolute) or  2 (relative)\n',...
         'noise.distr: 1 (uniform)  or  2 (Gauss)\n'])
  
noise = struct('noisefun',1,'level',0.01,'type',1,'distr',1)

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem definition

% define problem parameters (to be adapted to your problem)
n=5;  % dimension
p=2;  % Norm in objective function
e=1;  % Exponent in objective function function

% create random matrix and right hand side for objective function
% (specific to the model problem; replace by whatever data you
% need to provide to your objectiv function)
A=rand(n)-0.5; 
b=-sum(A,2);

z0 = ones(n,1);

% continuous bounds
clow = z0-10*ones(n,1); cupp = z0+10*ones(n,1); 

% create objective function f(x)=||Ax-b||_p^e
func  = @(z) norm(A*z-b,p).^e; 

% fun = @getf computes function value of func at z, while
% collecting statistics and enforcing stopping test
fun  = @getf;  

% To solve your own problem, simply replace in the previous line 
% the expression after @(z) by your expression or function call. 
% Parameters in this expression are passed by value, hence are 
% fixed during minimization.

% integer bounds
low = zeros(n,1); upp = 100*ones(n,1);
% starting point
z      = 50*ones(n,1);  

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tuning=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem with IMATRS

% pass stop and print criteria
% (indefinite run if no stopping criterion is given)
st = struct('secmax',180,'nfmax',500*n,'accf',0.001,...
           'fbest',0.01,'prt',1);

st.noise = noise; st.noisefun = noise.noisefun; st.z=z;

% initilize feinfo
initf(func,low,upp,clow,cupp,st)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  solve the problem with IMATRS
if Tuning % self-tuning and info
    % given are the defaults
    % only the deviating parameters need to be set!
    if feinfo.prt>0
      fprintf(dline)
      fprintf(['tune:   structure containing all tuning parameters\n',...
        dline,...
        'parameter for trust region condition: ',...
        'zeta = 1e-20 \n\n',...
        'parameter to update step size in line search: ',...
        'nu = 2 \n\n',...
        'parameter to update radius in trust region: ',...
        'theta= 2\n \n',...
        'initial trust region radius: ',...
        'Deltainit = 10 \n\n',...
        'pretty upper bound for the radius: ',...
        'Deltabar = 3 \n\n',...
        'mamimum value for sigma: sigmamax = 1e10 \n\n',...
        'factor for adjusting X: gammaX = 100\n \n',...
        'factor for adjusting sc and y: ',...
        'gammav = 100 \n',...
        dline]); 
    end
        tune = struct('zeta',1e-20,'nu',2,'theta',2,'Deltainit',10,...
            'sigmamax',1e10,'Deltabar',3,'gammaX',100,'gammav',100);
else
       tune = []; % full tuning inside IMATRS is used
end

% call solver
try
   tic
   IMATRS(fun,z,tune);
 catch ME
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % enforce stopping test %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(strfind(feinfo.error, 'not allowed'))
          % pass on error
          error(feinfo.error);
        end
        % if a solver fails
        % (which often happens because of a designed error in fun)
        if isempty(strfind(feinfo.error,'reached'))
          feinfo.error = ME.message;
          disp(' ');
          disp(['******************************************',...
                '*****************']);
          disp(['******************************************',...
                '*****************']);
          disp(['******************************************',...
                '*****************']);
          disp(['*** error: ',ME.message]);
          stack = ME.stack;
          for i = 1:length(stack)
            file = stack(i).file;
            name = stack(i).name;
            line = stack(i).line;
            disp(['*** line ',num2str(line),' in    ',name,...
                  '    from ',file]);
          end
          disp(['******************************************',...
                '*****************']);
          disp(['******************************************',...
                '*****************']);
          disp(['******************************************',...
                '*****************']);
          disp(' ');
        end
        disp(['stopped since ',feinfo.error]);
end
if feinfo.prt>0
  z  = feinfo.zbest     % point with best acc
end
f  = feinfo.fbest;     % function value at best acc
qf = feinfo.qfbest;   
nf = feinfo.nf;
maxInfeas = feinfo.maxInfeas;


nfree = sum(z > feinfo.low & z < feinfo.upp);
text1 = sprintf('f=%5.6e\n',f);
text2 = ['nf=',num2str(nf),', nfree=',num2str(nfree),...
        sprintf(', qf=%7.1e',qf')];
if maxInfeas > 0
   text3 = ['\n   maximal infeasibility of evaluation points: ',...
          num2str(maxInfeas)];
   text = [text1,text2,text3];
else
   text = [text1,text2];
end
text      = [text,', ', feinfo.error,'\n'];
showtime0 = showtime;
text      = [text,sprintf(showtime0),'\n','IMATRS completed - ',time];
text      = [dline,text,'\n',dline,'\n\n\n'];
% print on screen
if feinfo.prt>0, fprintf(1,text); end

for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
