r = 0.001:0.001:10.0 ;
l=9

%as =[0.2; 0.4; 0.6; 0.8; 1.0]% 4; 5; 6; 7; 8; 9; 10];
%as =[1.25;2.5;3.75; 5];
%as =[1;2;3;4;5];
as =[1;3;5];
%as =[0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0];
%as = [0.25; 0.5; 0.75; 1.0; 1.25; 1.5; 1.75; 2.0; 2.25; 2.5; 2.75; 3.0; 3.25; 3.5; 3.75; 4.0;4.25; 4.5;4.75; 5.0];

alphasFinal = [];
Betas = [];
Ys = [];
ngss = [];

for l=0:9
gs = [];
ngs = [];
alphas = [];


size(as)(1)
for i=1:size(as)(1)
check=0
while(check==0)
  alpha = (l+1)*20/i**1*rand();
  g = r.^l.*exp(-alpha*r.^2);
  ng = g/sqrt(trapz(r,r.*r.*g.*g));
  x = r(find(1 - cumtrapz(r,r.*r.*ng.*ng) < 0.01))(1);
%  [val, indx] = max(abs(r.*ng));
%  x = ng(indx);
  if( as(i) + 0.01  > x && x > as(i) - 0.01 )
    check=1;
    break;
end
end
  gs = [gs; g];
  ngs = [ngs; ng];
  alphas = [alphas alpha]
  x
end

%plot(r,ngs(1,:));

X = zeros(size(as)(1), size(as)(1));

for i=1:size(as)(1)
  for j=1:size(as)(1)
%    X(i,j) = trapz(r, r.*r.*ngs(i,:).*ngs(j,:))
    X(i,j) = trapz(r, r.*r.*gs(i,:).*gs(j,:))
  end
end

Beta = sqrtm(inv(X));

%Y = Beta*ngs;
Y = Beta*gs;

for i=1:size(as)(1)
  for j=1:size(as)(1)
    trapz(r,r.*r.*Y(i,:).*Y(j,:))
  end
end

for i=1:size(as)(1)
plot(r,r.*r.*Y(1,:).*Y(1,:),'r')
axis([0 10 -2.5 4])
hold on
plot(r,r.*r.*Y(2,:).*Y(2,:),'b')
plot(r,r.*r.*Y(3,:).*Y(3,:),'g')
%plot(r,r.*r.*Y(4,:).*Y(4,:),'r')
%plot(r,r.*r.*Y(5,:).*Y(5,:),'b')
end
Betas = cat(3,Betas,Beta);
Ys = cat(3,Ys,Y);
ngss = cat(3,ngss,ngs);
alphasFinal = [alphasFinal; alphas];
%printf("HALLELOYA ")
l
end

save parameters3.dat *
print -dpsc
print -deps -color fo.eps


%#plot(r,r.*r.*ng3.*ng3);
%pause(10);
