#!/usr/bin/octave
# vim:ft=matlab

% sdp-soln
values
n = size(A,1);

% omega
O = 2*eye(n); O(1) = 1; O(end) = 1;
Op = O^0.5;
Om = O^-0.5;


org = Op * A * Om;
sym = Op * A * Om;

for i=1:size(A,1)
	for j=i:size(A,2)
		sym(j,i) = sym(i,j);
	end
end

C = sym;

res   = [];
yvals = [];
rv    = [];
xs    = [];
ys    = [];
for i=1:1000
	r = 2*rand(size(A,1),1)-1;
	r = r / norm(r);

	% project v_i onto line defined by r
	y = V'*r;
	%y = y / norm(y);
	y = sort(y);
	ys = [ ys y ];

	x = Op \ y;
	yvals = [ yvals y'*y ];
	xs    = [ xs x ];

	val = y'*C*y;
	res = [ res val ];

	rv  = [ rv r'*C*r ];
	R = r*r';
end

v = V(1,:)';

[max_y2,i]=max(yvals)
[max_obj,i]=min(res)
max_random=min(rv)
v_1_obj=v'*C*v

x = xs(:,i);
y = ys(:,i);

done = zeros(n,1);
path = [];

x = x + min(x) + 1;
[m,i] = max(x);
while(size(path) < n)
	done(i) = 1;
	path = [path i];
	x(i) = -1;
	j    = i;
	[m,i] = max(x.*A(:,i));
	if( (done(i) || A(i,j) <1 ) && size(path) < n)
		[m,i] = max(x);
		printf('Jumping from %d to %d...\n', j,i);
	end
end
