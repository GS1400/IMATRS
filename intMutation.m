%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% intMutation.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform integer mutation

function [xbest,fbest,xf,a,good,allones,succ,OffspringPop,restart]=...
         intMutation(fun,n,xbest,fbest,a,succ,allones,D,nx,M,xf,tune)

global feinfo

good=0; OffspringPop=[]; restart=0;
for l=1:size(D,2)
    Offspring.sd = D(:,l);
    Offspring.p  = round(M*Offspring.sd);
    p =Offspring.p;
    if norm(p)==0, p = D(:,l); restart=1; end
    % initialize vector alpha_max
   alphamax = Inf * ones(n,1);
    % caluclate max alpha
   ind = (p > 0 );
   alphamax(ind)=( feinfo.upp(ind) - xbest(ind) )./ p(ind);
   ind = ( p < 0 );
   alphamax(ind)=( feinfo.low(ind) - xbest(ind) )./ p(ind);
   % compute starting alpha
   alpha_bar  = round( min(alphamax) );
   if a(l)<alpha_bar
      alphain    = round(min(a(l),sqrt(a(l)*alpha_bar)));
   else
      alphain    = round(min(a(l),alpha_bar));
   end
  [alpha,ytrial,ftrial,xf]=intLSS(fun,nx,xf,xbest,...
                     p,alphain,alpha_bar,fbest,tune);
  if feinfo.done; return; end 

  if alpha<=0
    p=-p; 
    % initialize vector alpha_max
    alphamax = Inf * ones(n,1);
    % caluclate max alpha
    ind = ( p > 0 );
    alphamax(ind)=( feinfo.upp(ind) - xbest(ind) )./ p(ind);
    ind = ( p < 0 );
    alphamax(ind)=( feinfo.low(ind) - xbest(ind) )./ p(ind);
    %compute starting alpha
    alpha_bar  = round( min(alphamax));

    if a(l)<alpha_bar
        alphain    = round(min(a(l),sqrt(a(l)*alpha_bar)));
    else
         alphain    = round(min(a(l),alpha_bar));
    end

    [alpha,ytrial,ftrial,xf]=intLSS(fun,nx,xf,xbest,...
     p,alphain,alpha_bar,fbest,tune);

     if feinfo.done; return; end 

     if alpha>0
          Offspring.sd = -Offspring.sd;
          Offspring.p = - Offspring.p;
          succ(l) = succ(l)+1;
          if allones >= 1
              allones = 0;
          end
         a(l) = alpha;
          Offspring.f = ftrial; 
          OffspringPop= [OffspringPop Offspring];
          if(ftrial < fbest)
             fbest = ftrial;
             xbest = ytrial; good=1;
             if feinfo.prt>=3
               fprintf('integer line search was successful\n')
             end
          end
     else
          a(l) = max(1,floor(a(l)/tune.nu));
     end
  else
      succ(l) = succ(l)+1;
      if allones >= 1, allones = 0; end
      a(l) = alpha;
      if(ftrial < fbest)
          fbest = ftrial; xbest = ytrial; good=1;
          if feinfo.prt>=3
             fprintf('integer line search was successful\n')
          end
      end
      Offspring.f = ftrial;
      OffspringPop= [OffspringPop Offspring];
  end

  if (allones > 1), break; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%