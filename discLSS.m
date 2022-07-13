%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% discLSS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform integer line search
%
% discLSS is a special case of nm_discrete_linesearch from box_DFL 
% by G.Liuzzi, S.Lucidi, F.Rinaldi

function [alpha,x,f,xf] = discLSS(fun,nx,xf,y,d,...
                                          alphain,alpha_bar,fb,tune)
global feinfo
x=y; f=fb; alpha = [];
% calculate dimension of the problem
n = length(d);
% Build first point for starting linesearch
if (alphain > 0)
    ytrial = y + alphain * d;
    sxf  = size(xf,1);
    diff = (xf(:,1:n)-repmat(ytrial',sxf,1)).^2;
    [mn,ind]=min(sum(diff,2));
    if (mn<=10^-16)
        fval   = xf(ind(1),n+1); ftrial = fval;
    else
        [ftrial]= fun(ytrial);
        xf(feinfo.nf,1:nx)=ytrial'; xf(feinfo.nf,nx+1)=ftrial;
        if feinfo.done; return; end 
    end
else
    ftrial = Inf;
end
% cicle for updating alpha
if (alphain > 0) && (ftrial<fb)
    % initialize alpha and best point
    alpha=alphain; x=ytrial; f=ftrial;
    % calculate trial point
    if(alpha < alpha_bar)
        ytrial = y + min(alpha_bar,tune.nu*alpha)*d;
        sxf=size(xf,1);
        diff=(xf(:,1:n)-repmat(ytrial',sxf,1)).^2;
        [mn,ind]=min(sum(diff,2));
        if (mn<=10^-16)
            fval = xf(ind(1),n+1); ftrial = fval;
        else
            [ftrial] = fun(ytrial);
            xf(feinfo.nf,1:nx)=ytrial';
            xf(feinfo.nf,nx+1)=ftrial;
            if feinfo.done; return; end 
        end
    else
        ftrial = Inf;
    end
    % expansion step (increase stepsize)
    while (alpha<alpha_bar) && (ftrial < fb)
        % alpha calulation and best point updating
        alpha=min(alpha_bar, tune.nu*alpha);
        % best point updating
        x=ytrial; f=ftrial;
        % next point to be tested
        if(alpha < alpha_bar)
            ytrial = y + min(alpha_bar,tune.nu*alpha)*d;
            sxf=size(xf,1);
            diff=(xf(:,1:n)-repmat(ytrial',sxf,1)).^2;
            [mn,ind]=min(sum(diff,2));
            if (mn<=10^-16)
               fval = xf(ind(1),n+1); ftrial = fval;
            else
               [ftrial] = fun(ytrial);
               xf(feinfo.nf,1:nx)=ytrial'; xf(feinfo.nf,nx+1)=ftrial;
               if feinfo.done; return; end 
            end
        else
            ftrial = Inf;
        end
    end
else
    alpha = 0; x = y; f = Inf;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%