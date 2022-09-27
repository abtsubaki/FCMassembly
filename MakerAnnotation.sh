#!/bin/bash
#PBS -N MAKER
#PBS -e maker.err
#PBS -o maker.out
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=745:00:00
#PBS -l mem=100gb
#PBS -m abe

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8


module load app/maker/maker

#provide maker with annotation details: see maker_opts.log for options
echo "creating ctl files\n"
maker -CTL
echo "completed ctl files\n"

echo "running maker\n"
maker
echo "completed maker\n"

#merge gff output files
echo "running merge gff\n"
gff3_merge -s -d $DIR/fcmgenomev2_purge.maker.output/fcmgenomev2_purge_master_datastore_index.log > $DIR/fcmgenome_purge_maker_mix.gff
echo "done running merge gff\n"

