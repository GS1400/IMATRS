%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% prime_vector.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether d is a prime vector
%
% prime_vector comes from box_DFL by G.Liuzzi, S.Lucidi, F.Rinaldi

function flag = prime_vector(d)
n = length(d);
flag = 0;
if(n==1), flag = 1;return; end
temp = gcd(abs(d(1)),abs(d(2)));
if(n==2), flag = (temp == 1); return; end
for i = 3:n
    temp = gcd(temp,abs(d(i)));
    if temp == 1, flag = 1;return; end
end
if temp ~= 1, flag = 0; return; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 