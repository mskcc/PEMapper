#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/bin:$PATH
source $SDIR/bin/lsf.sh

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
PIPENAME="PEMapper"

##
# Process command args

TAG=qPEMAP

COMMAND_LINE=$*
function usage {
    echo
    echo "usage: $PIPENAME/pipe.sh [-s SAMPLENAME] GENOME SAMPLEDIR"
    echo "version=$SCRIPT_VERSION"
    echo "    -g ListGenomes"
    echo
    exit
}

BWA_OPTS=""
SAMPLENAME="__NotDefined"
while getopts "s:hgb:" opt; do
    case $opt in
        s)
            SAMPLENAME=$OPTARG
            ;;
        b)
            BWA_OPTS=$BWA_OPTS" -"$OPTARG
            ;;
        h)
            usage
            ;;
        g)
            echo Currently defined genomes
            echo
            ls -1 $SDIR/lib/genomes
            echo
            exit
            ;;
        \?)
            usage
            ;;
    esac
done

shift $((OPTIND - 1))
if [ "$#" -lt "2" ]; then
    usage
fi

BWA_OPTS=$(echo $BWA_OPTS | perl -pe "s/^\s+//")
echo BWA_OPTS="["$BWA_OPTS"]"

GENOME=$1
shift

if [ -e $SDIR/lib/genomes/$GENOME ]; then
    source $SDIR/lib/genomes/$GENOME
else
    if [ -e $GENOME ]; then
        source $GENOME
    else
        echo
        echo GENOME=$GENOME Not Defined
        echo "Currently available (builtin) genomes"
        ls -1 $SDIR/lib/genomes
        echo
        exit
    fi
fi

SAMPLEDIR=$1
SAMPLEDIR=$(echo $SAMPLEDIR | sed 's/\/$//' | sed 's/;.*//')

SAMPLEDIRS=$*
SAMPLEDIRS=$(echo $SAMPLEDIRS | tr ';' ' ')

if [ $SAMPLENAME == "__NotDefined" ]; then
    SAMPLENAME=$(basename $SAMPLEDIR)
    if [ "$SAMPLENAME" == "" ]; then
        echo "Error in sample name processing; Null sample name"
        exit
    fi
fi

echo SAMPLENAME=$SAMPLENAME
TAG=${TAG}_$$_$SAMPLENAME

export SCRATCH=$(pwd)/_scratch/$(uuidgen -t)
mkdir -p $SCRATCH
echo SAMPLENAME=$SAMPLENAME >> $SCRATCH/RUNLOG
echo BWA_OPTS=$BWA_OPTS >> $SCRATCH/RUNLOG
echo GENOME=$GENOME >> $SCRATCH/RUNLOG
echo TAG=$TAG >> $SCRATCH/RUNLOG


##
# HiSeq TrueSeq maximal common adapter

ADAPTER="AGATCGGAAGAGC"
BWA_VERSION=$(bwa 2>&1 | fgrep Version | awk '{print $2}')

JOBS=""
BAMFILES=""

FASTQFILES=$(find -L $SAMPLEDIRS -name "*[_.]R1*.fastq.gz")
echo "FASTQFILES="$FASTQFILES

if [ "$FASTQFILES" == "" ]; then
    echo "Can not find any FASTQFILES"
    exit
fi


for FASTQ1 in $FASTQFILES; do

	case "$FASTQ1" in
		*_R1_*)
		FASTQ2=${FASTQ1/_R1_/_R2_}
		;;

		*.R1.*)
		FASTQ2=${FASTQ1/.R1./.R2.}
		;;

		*)
		echo
		echo "FATAL ERROR; INVALID FASTQ1 filename =" $FASTQ1
		exit

	esac

    BASE1=$(echo $FASTQ1 | tr '/' '_')
    BASE2=$(echo $FASTQ2 | tr '/' '_')
    UUID=$(uuidgen)

    # if MINLENGTH not set in ENV then set to 1/2 read length
    if [ "$MINLENGTH" == "" ]; then

        # Get readlength
        ONE_HALF_READLENGTH=$(zcat $FASTQ1 | $SDIR/bin/getReadLength.py | awk '{printf("%d\n",$1/2)}')
        echo ONE_HALF_READLENGTH=$ONE_HALF_READLENGTH
        echo ONE_HALF_READLENGTH=$ONE_HALF_READLENGTH >> $SCRATCH/RUNLOG
        export MINLENGTH=$ONE_HALF_READLENGTH

    fi

    QRUN 2 ${TAG}_MAP_01__$UUID VMEM 5 \
        clipAdapters.sh $ADAPTER $FASTQ1 $FASTQ2
    CLIPSEQ1=$SCRATCH/${BASE1}___CLIP.fastq
    CLIPSEQ2=$SCRATCH/${BASE2}___CLIP.fastq

    BWA_THREADS=8

    echo -e "@PG\tID:$PIPENAME\tVN:$SCRIPT_VERSION\tCL:$0 ${COMMAND_LINE}" >> $SCRATCH/${BASE1%%.fastq*}.sam

    QRUN $BWA_THREADS ${TAG}_MAP_02__$UUID HOLD ${TAG}_MAP_01__$UUID VMEM 8 \
        bwa mem $BWA_OPTS -t $BWA_THREADS $GENOME_BWA $CLIPSEQ1 $CLIPSEQ2 \>\>$SCRATCH/${BASE1%%.fastq*}.sam

    QRUN 4 ${TAG}_MAP_03__$UUID HOLD ${TAG}_MAP_02__$UUID VMEM 33 \
        picardV2 AddOrReplaceReadGroups MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true SO=coordinate \
        LB=$SAMPLENAME PU=${BASE1%%_R1_*} SM=$SAMPLENAME PL=illumina CN=GCL \
        I=$SCRATCH/${BASE1%%.fastq*}.sam O=$SCRATCH/${BASE1%%.fastq*}.bam

    BAMFILES="$BAMFILES $SCRATCH/${BASE1%%.fastq*}.bam"

done

echo
echo BAMFILES=$BAMFILES
echo HOLDTAG="${TAG}_MAP_*"
echo BAMFILES=$BAMFILES >> $SCRATCH/RUNLOG
echo HOLDTAG="${TAG}_MAP_*" >> $SCRATCH/RUNLOG
echo

INPUTS=$(echo $BAMFILES | tr ' ' '\n' | awk '{print "I="$1}')

BWATAG=$(echo $BWA_OPTS | perl -pe 's/-//g' | tr ' ' '_')

OUTDIR=out___$BWATAG
mkdir -p $OUTDIR
QRUN 4 ${TAG}__04__MERGE HOLD "${TAG}_MAP_*"  VMEM 33 LONG \
    picardV2 MergeSamFiles SO=coordinate CREATE_INDEX=true \
    O=$OUTDIR/${SAMPLENAME}.bam $INPUTS
    
QRUN 4 ${TAG}__05__MD HOLD ${TAG}__04__MERGE VMEM 33 LONG \
    picardV2 MarkDuplicates \
    I=$OUTDIR/${SAMPLENAME}.bam \
    O=$OUTDIR/${SAMPLENAME}___MD.bam \
    M=$OUTDIR/${SAMPLENAME}___MD.txt \
    CREATE_INDEX=true \
    R=$GENOME_FASTA
    
QRUN 4 ${TAG}__05__STATS HOLD ${TAG}__05__MD VMEM 33 LONG \
    picardV2 CollectAlignmentSummaryMetrics \
    I=$OUTDIR/${SAMPLENAME}___MD.bam O=$OUTDIR/${SAMPLENAME}___AM.txt \
    R=$GENOME_FASTA \
    LEVEL=null LEVEL=SAMPLE

QRUN 4 ${TAG}__05__STATS HOLD ${TAG}__05__MD VMEM 33 LONG \
    picardV2 CollectInsertSizeMetrics \
    I=$OUTDIR/${SAMPLENAME}___MD.bam O=$OUTDIR/${SAMPLENAME}___INS.txt \
	H=$OUTDIR/${SAMPLENAME}___INSHist.pdf \
    R=$GENOME_FASTA

QRUN 4 ${TAG}__05__STATS HOLD ${TAG}__05__MD VMEM 33 LONG \
    picardV2 CollectGcBiasMetrics \
    I=$OUTDIR/${SAMPLENAME}___MD.bam O=$OUTDIR/${SAMPLENAME}___GCB.txt \
    CHART=$OUTDIR/${SAMPLENAME}___GCB.pdf \
    S=$OUTDIR/${SAMPLENAME}___GCBsummary.txt \
    R=$GENOME_FASTA

QRUN 4 ${TAG}__05__STATS HOLD ${TAG}__05__MD VMEM 33 LONG \
    picardV2 CollectWgsMetrics \
    I=$OUTDIR/${SAMPLENAME}___MD.bam O=$OUTDIR/${SAMPLENAME}___WGS.txt \
    R=$GENOME_FASTA

#QRUN 4 ${TAG}__05__DOWN HOLD ${TAG}__04__MERGE VMEM 33 LONG \
#    picardV2 DownsampleSam \
#    I=$OUTDIR/${SAMPLENAME}.bam \
#    O=$OUTDIR/${SAMPLENAME}___Dn10.bam \
#    P=0.1 CREATE_INDEX=true

if [ "$DBSNP" != "" ]; then
    QRUN 4 ${TAG}__05__OXO HOLD ${TAG}__05__MD VMEM 33 LONG \
        picardV2  CollectOxoGMetrics \
        R=$GENOME_FASTA \
        DB_SNP=$DBSNP \
        I=$OUTDIR/${SAMPLENAME}___MD.bam \
        O=$OUTDIR/${SAMPLENAME}___OxoG.txt
else
    QRUN 4 ${TAG}__05__OXO HOLD ${TAG}__05__MD VMEM 33 LONG \
        picardV2  CollectOxoGMetrics \
        R=$GENOME_FASTA \
        I=$OUTDIR/${SAMPLENAME}___MD.bam \
        O=$OUTDIR/${SAMPLENAME}___OxoG.txt
fi

QRUN 1 ${TAG}__07_CLEANUP HOLD ${TAG}__05__STATS \
    rm -rf $SCRATCH

