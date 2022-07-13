%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intTRS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs an integer trust region algorithm

function [xbest,fbest,good,xf]=intTRS(fun,n,xbest,fbest,nx,xf,M,tune)
global feinfo
Deltadis  = tune.Deltainit; good=0; ii=0; jj=0;
while 1
    ii = ii+1;
    X = xf(feinfo.nf-n:feinfo.nf,1:nx)';
    F = xf(feinfo.nf-n:feinfo.nf,nx+1)';
    gz = fitGrad(n,X,F,xbest,fbest,n+1,tune);
    M = adjustY(M,tune); 
    warning off
    invM = inv(M);
    warning on
    ll = max(feinfo.low-xbest,-Deltadis);
    uu = min(feinfo.upp-xbest,Deltadis);
    igls = 0;
    for i = 1 : n
        if ll(i) >= uu(i), igls=1; end
    end
    if ~igls && rank(M) == n
       r = -invM*gz; sz = sils(M,r,1);
       sz =max(ll,min(sz,uu));
       if isempty(sz)||norm(sz)==0, break; end
    else
        break
    end
    ytrial = min(feinfo.upp,max(feinfo.low,xbest+sz));
    sxf=size(xf,1);
    diff=(xf(:,1:n)-repmat(ytrial',sxf,1)).^2;
    [mn,~]=min(sum(diff,2));
    if (mn<=10^-16)
        jj=jj+1; Deltadis=Deltadis-2;
        if Deltadis<=0, break; end
        if jj==tune.Deltainit, break ;end
        continue
    else
       [ftrial] = fun(ytrial);
       xf(feinfo.nf,1:nx)=ytrial';
       xf(feinfo.nf,nx+1)=ftrial;
       if feinfo.done; return; end 
    end          
    gzsz = max(-gz'*sz,eps*abs(gz)'*abs(sz));
    mueff = (fbest-ftrial)/gzsz;
    succTR  = (mueff*abs(mueff-1) > tune.zeta);
    % Updating iterate and trust-region radius.
    if succTR
       if feinfo.prt>=3
          fprintf('integer trust region was successful\n')
       end 
       ytrial0  = ytrial; ftrial0 = ftrial; d=sz;
       % initialize vector alpha_max
       alphamax = Inf * ones(n,1);
       % caluclate max alpha
       ind = ( d > 0 );
       alphamax(ind)=( feinfo.upp(ind) - xbest(ind) )./ d(ind);
       ind = ( d < 0 );
       alphamax(ind)=( feinfo.low(ind) - xbest(ind) )./ d(ind);
       % compute starting alpha
       alpha_bar  = round( min(alphamax));
       alphain    = min(1,alpha_bar);
       [alpha,ytrial,ftrial,xf]= ...
                   intLSS(fun,nx,xf,xbest,...
                   d,alphain,alpha_bar,fbest,tune);
       if feinfo.done; return; end 
       if alpha<0
            d=-d;
            % initialize vector alpha_max
            alphamax = Inf * ones(n,1);
            % caluclate max alpha
           ind = ( d > 0 );
           alphamax(ind)=( feinfo.upp(ind) - xbest(ind) )./ d(ind);
           ind = ( d < 0 );
           alphamax(ind)=( feinfo.low(ind) - xbest(ind) )./ d(ind);
           % compute starting alpha
           alpha_bar  = round( min(alphamax));
           alphain   = min(1,alpha_bar);
           [alpha,ytrial,ftrial,xf]=...
                          intLSS(fun,nx,xf,xbest,p,...
                          alphain,alpha_bar,fbest,tune);
          if feinfo.done; return; end                
          if alpha>0
             if(ftrial < fbest)
               fbest = ftrial; xbest = ytrial; 
               if feinfo.prt>=3
                  fprintf('integer line search was successful\n')
               end
             else
                xbest  = ytrial0; fbest = ftrial0;
             end
          else
             xbest  = ytrial0; fbest = ftrial0; 
          end
      else
           if alpha>0
             if(ftrial < fbest)
                if feinfo.prt>=3
                   fprintf('integer line search was successful\n')
                end
                fbest = ftrial; xbest = ytrial; 
             else
                xbest  = ytrial0; fbest = ftrial0; 
             end
           else
               xbest  = ytrial0; fbest = ftrial0;  
           end
      end
     good=1;
     break
    else
        good=0;
        if Deltadis<=tune.Deltabar,Deltadis=Deltadis-1;
        else, Deltadis = floor(Deltadis/tune.theta);
        end
        if feinfo.prt>=3
           fprintf('integer trust region was unsuccessful\n')
       end 
    end
    if Deltadis<=0, break; end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%