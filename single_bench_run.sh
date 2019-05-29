BENCHFOLDER="./benchmarks"
defaultActivity=0.5
nvec=10000
ninst=1000
theta=0.1
N=1000
trig=8

echo "./ga $BENCHFOLDER/$1.v $nvec $theta $trig $ninst $1 $2 $N 0 0 $defaultActivity -MERO"
./ga $BENCHFOLDER/$1.v $nvec $theta $trig $ninst $1 $2 $N 0 0 $defaultActivity -MERO

sleep 10

echo "./ga $BENCHFOLDER/$1.v $nvec $theta $trig $ninst $1 $2 $N 0 0 $defaultActivity -genePatternPairGA"
./ga $BENCHFOLDER/$1.v $nvec $theta $trig $ninst $1 $2 $N 0 0 $defaultActivity -genePatternPairGA
