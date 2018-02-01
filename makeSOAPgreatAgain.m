r=0.00005:0.00005:20.0;
l=9

as =[1; 2; 3; 4; 5]% 4; 5; 6; 7; 8; 9; 10];

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
  alpha = 20*rand();
  g = r.^l.*exp(-alpha*r.^2);
  ng = g/sqrt(trapz(r,r.*r.*g.*g));
  x = r(find(1 - cumtrapz(r,r.*r.*ng.*ng) < 0.001))(1) ;
  if( as(i)+0.001  > x && x > as(i)-0.001 )
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

X = zeros(size(as)(1), size(as)(1))

for i=1:size(as)(1)
  for j=1:size(as)(1)
    X(i,j) = trapz(r, r.*r.*ngs(i,:).*ngs(j,:))
  end
end

Beta = sqrtm(inv(X))

Y = Beta*ngs;

for i=1:size(as)(1)
  for j=1:size(as)(1)
    trapz(r,r.*r.*Y(i,:).*Y(j,:))
  end
end

for i=1:size(as)(1)
plot(r(1:end/2),r(1:end/2).*r(1:end/2).*Y(1,:)(1:end/2).*Y(1,:)(1:end/2),'r')
axis([0 10 -2.5 4])
hold on
plot(r(1:end/2),r(1:end/2).*r(1:end/2).*Y(2,:)(1:end/2).*Y(2,:)(1:end/2),'b')
plot(r(1:end/2),r(1:end/2).*r(1:end/2).*Y(3,:)(1:end/2).*Y(3,:)(1:end/2),'g')
plot(r(1:end/2),r(1:end/2).*r(1:end/2).*Y(4,:)(1:end/2).*Y(4,:)(1:end/2),'r')
plot(r(1:end/2),r(1:end/2).*r(1:end/2).*Y(5,:)(1:end/2).*Y(5,:)(1:end/2),'b')
end
Betas = cat(3,Betas,Beta);
Ys = cat(3,Ys,Y);
ngss = cat(3,ngss,ngs);
alphasFinal = [alphasFinal; alphas];
%printf("HALLELOYA ")
l
end

save parameters.dat *
print -dpsc
print -deps -color fo.eps


%#plot(r,r.*r.*ng3.*ng3);
%pause(10);
