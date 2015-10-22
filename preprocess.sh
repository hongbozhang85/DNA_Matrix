#! /bin/bash

# if [feedback != 1] then no debug information on terminal
feedback=1
#cluster=2
#motiflen=9

# set paths and filenames
fastaname=$1
ininame=$fastaname.ini
datname=$fastaname.dat

cd data/
datapath=$PWD
cd ../input/
inputpath=$PWD
cd ../output/
outputpath=$PWD
cd ../

fastafile=$datapath/$fastaname
inifile=$inputpath/$ininame
datfile=$inputpath/$datname

if [ $feedback -eq 1 ]
then
echo 'data path: ' $datapath
echo 'input path: ' $inputpath
echo 'output path: ' $outputpath
echo 'current fasta file: ' $fastafile
echo 'ini file: ' $inifile
echo 'dat file: ' $datfile
fi

# generate *.dat in /input/
if [ -e $datfile ]
then
rm $datfile
fi
sed '1d' $fastafile | tr -d "[N\n]" > $datfile
echo $datname 'has been generated'

#generate *.ini in /input/
if [ -e $inifile ]
then
rm $inifile
fi
sed -n '1p' $fastafile > $inifile
wc -m $datfile >> $inifile
echo $datfile >> $inifile
echo $inputpath >> $inifile
echo $datname >> $inifile
echo $outputpath >> $inifile
#echo $cluster >> $inifile
#echo $motiflen >> $inifile
echo $ininame 'has been generated'

#run ./genomics xxx.ini
./genomics $inifile

##run ./image xxx.ini
#./image $inifile

##run ./motscurve xxx.ini
#./motscurve $inifile

##run ./strscurve xxx.ini
#./strscurve $inifile
