function B = rotateXYZ(A,myalpha, mybeta, mygamma)

  Rx = zeros(3,3);
  Ry = zeros(3,3);
  Rz = zeros(3,3);
  rotationMat = zeros(3,3) ;

  B = zeros(size(A)(1),size(A)(2));

  Rx(1,1)= 1; Rx(1,2)= 0;            Rx(1,3)= 0;
  Rx(2,1)= 0; Rx(2,2)= cos(myalpha); Rx(2,3)= -sin(myalpha);
  Rx(3,1)= 0; Rx(3,2)= sin(myalpha); Rx(3,3)= cos(myalpha);

  Ry(1,1) = cos(mybeta);   Ry(1,2) = 0;  Ry(1,3) = sin(mybeta);
  Ry(2,1) = 0;             Ry(2,2) = 1;  Ry(2,3) = 0;
  Ry(3,1) = -sin(mybeta);  Ry(3,2) = 0;  Ry(3,3) = cos(mybeta);

  Rz(1,1) = cos(mygamma);  Rz(1,2) = -sin(mygamma); Rz(1,3) = 0;
  Rz(2,1) = sin(mygamma);  Rz(2,2) = cos(mygamma);  Rz(2,3) = 0;
  Rz(3,1) = 0;             Rz(3,2) = 0;             Rz(3,3) = 1;

  rotationMat = Rx*Ry*Rz;


  for i=1:size(A)(1) 
      B(i,:) = A(i,:)*rotationMat; 
  end

end
