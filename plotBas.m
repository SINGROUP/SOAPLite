
alphas = load("alphasPy.dat")
betas = load("betasPy.dat") %only for l=0
r = [0.0:0.01:10];

g1 = zeros(1,1001);
g2 = zeros(1,1001);
g3 = zeros(1,1001);
g4 = zeros(1,1001);
g5 = zeros(1,1001);
g6 = zeros(1,1001);
g7 = zeros(1,1001);
g8 = zeros(1,1001);
g9 = zeros(1,1001);
g10 = zeros(1,1001);
betas = reshape(betas, [10,10])

for k=1:10
  g1 = g1 + betas(1,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g2 = g2 + betas(2,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g3 = g3 + betas(3,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g4 = g4 + betas(4,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g5 = g5 + betas(5,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g6 = g6 + betas(6,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g7 = g7 + betas(7,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g8 = g8 + betas(8,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g9 = g9 + betas(9,k)*exp(-alphas(k)*r.*r);
end
for k=1:10
  g10 = g10 + betas(10,k)*exp(-alphas(k)*r.*r);
end

plot(r,g1)
hold on
plot(r,g2,'r')
plot(r,g3,'g')
plot(r,g4,'b')
plot(r,g5,'b')
plot(r,g6,'r')
plot(r,g7,'g')
plot(r,g8,'r')
plot(r,g9,'b')
plot(r,g10,'r')
xlabel("r[A]")


print -dpsc
print -deps -color fo.eps


pause(5)
