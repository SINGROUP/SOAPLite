../Bash/getXYZ.sh $1 $2 > buff.xyz
../Bash/addH.sh buff.xyz;
a=$(wc -l < buff.xyz);
b=$(wc -l < $3);
c=$(($b - 1))
echo $a;
mkdir Dir_$2
for i in `seq 1 $a`; do mkdir Dir_$2/XYZ_$i; done
for i in `seq 1 $a`; do cat $3 > Dir_$2/XYZ_$i/$3.$i.xyz; done
for i in `seq 1 $a`; do sed -n "${i}p" buff.xyz >> Dir_$2/XYZ_$i/$3.$i.xyz; done
for i in `seq 1 $a`; do sed -i "s/^\s*[0-9][0-9]*\s*$/$c/g" Dir_$2/XYZ_$i/$3.$i.xyz ; done 
