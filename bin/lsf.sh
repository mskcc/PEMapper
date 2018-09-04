QSYNC=#


##
# QRUN ALLOC QTAG <HOLD hold_id> <VMEM size> LONG
#
# HOLD, VMEM and LONG are optional but if you more than one
# they must be in this order
#   HOLD
#   VMEM
#   LONG
#

QRUN () {

    ALLOC=$1
    QTAG=$2
    echo QTAG=$QTAG
    shift 2
    QHOLD=""
    if [ "$1" == "HOLD" ]; then
        QHOLD="-w post_done($2)"
        shift 2
        echo QHOLD=$QHOLD
    fi

    VMEM=""
    if [ "$1" == "VMEM" ]; then
        VMEM='-R "rusage[mem='$2']"'
        shift 2
        echo VMEM=$VMEM
    fi

    TIME="-We 119"
    if [ "$1" == "LONG" ]; then
        TIME=""
        shift 1
        echo LONG Job
    fi

    RET=$(bsub $TIME $QHOLD $VMEM -n $ALLOC -J $QTAG -o LSF.PEMAP/ $*)
    echo RET=bsub $QHOLD $VMEM -n $ALLOC -J $QTAG -o LSF.PEMAP/ $*
    echo "#QRUN RET=" $RET
    echo
    JOBID=$(echo $RET | perl -ne '/Job <(\d+)> /;print $1')

}
