n = 5
m = 5

x = [-6:0.02:6];
y = [-6:0.02:6];
z = [-6:0.02:6];

[X,Y,Z] = meshgrid(x,y,z);

a = load("typeCu.dat");
b = load("betas.dat");
alphas = load("alphas.dat");

b = reshape(b,[5 5 10]);

x0 = a(:,1);
y0 = a(:,2);
z0 = a(:,3);

x0 = x0 - 3.24644;
y0 = y0 + 2.62429;
z0 = z0 - 15.5668;

size(a)

G1 = b(n,1,1).*exp(-alphas(1,1).* (X .^2 + Y.^2 + Z.^2));
G2 = b(m,1,1).*exp(-alphas(1,1).* (X .^2 + Y.^2 + Z.^2));
for k =2:5
  G1 = G1 + b(n,k,1).*exp(-alphas(k,1).* ((X).^2 + (Y).^2 + (Z).^2));
  G2 = G2 + b(m,k,1).*exp(-alphas(k,1).* ((X).^2 + (Y).^2 + (Z).^2));
end
F1 =  G1.*exp(-(X - x0(1)).^2 - (Y - y0(1)).^2 - (Z - z0(1)).^2);
F2 =  G2.*exp(-(X - x0(1)).^2 - (Y - y0(1)).^2 - (Z - z0(1)).^2);
for i=2:size(a)(1)
  F1 = F1 + G1.*exp(-(X - x0(i)).^2 - (Y - y0(i)).^2 - (Z - z0(i)).^2);
  F2 = F2 + G2.*exp(-(X - x0(i)).^2 - (Y - y0(i)).^2 - (Z - z0(i)).^2);
end
final1=  trapz(z,trapz(y,trapz(x,F1)));
final2=  trapz(z,trapz(y,trapz(x,F2)));

final = final1*final2/4/pi
