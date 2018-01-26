a=$(wc -l < $1)

for i in `seq 1 $a`; do b=$(sed -n "${i}p" $1); sed -n "${b}p" $2; done
