#! /bin/bash                                                                   
#PBS -N matlab
#PBS -m a 
#PBS -m e
#PBS -M lengentyh@gmail.com  
##PBS -q workq
##---------------------------------------------------------------------------------------------------
## -l select=<# of chunks>[:<ncpus=#>], how many chunks. ncpus, how many cpus for each chunk.
## chunks is the set of resources allocated, typically, there is one chunk per MPI process.
## recommend select=2:ncpus=1 for 2-rank MPI job, select=1:ncpus=2 for 2-threads openMP job.
## select=<...>:cpu=6130 or cpu=2690v4 for request the 6130 CPU nodes, E5-2690v4 CPU nodes.
##---------------------------------------------------------------------------------------------------
#PBS -l select=1:ncpus=1
##---------------------------------------------------------------------------------------------------
## -l place= <type>[:<sharing>][:<group>]
## place the chunks on which way.
## vnodes is a virtual node, a abstract object representing a sub set of resources of a machine(host).
###     free,     on any vnode(s)
###     pack,     all chunks will be taken from one host
###     scatter,  only one chunk is taken from any host(means try to distribute your MPI rank to nodes)
###     vscatter, only one chunk is taken from any vnode(same as above but to vnodes)
## recommned place=pack for avoiding the hosts-crossing jobs
## place=pack:excl(exclhost) for exclusively using the vnode(node)
## a placement set will try to allocate the node for the jobs best,
## example: <...>:group=switch will arrange the vnodes(or nodes) for jobs at the same switch if possible
##---------------------------------------------------------------------------------------------------
#PBS -l place=pack:group=switch
##---------------------------------------------------------------------------------------------------

if ! $(type module &>/dev/null) || [[ -z "$LMOD_CMD" ]]; then
    echo "The \`module' or lua mod command failed, pleae contact the administrator"
    exit
fi

module load matlab

test -d "$PBS_O_WORKDIR" && cd "$PBS_O_WORKDIR" && echo -e "you are acquiring the resoures\n$(cat $PBS_NODEFILE)"

# user /tmp/$USER as the local scratch
export JOB_SCRATCH=/tmp/${USER}/${PBS_JOBID}
mkdir -p ${JOB_SCRATCH} 


use_record=F  # change to T if you want to recored the cpu and memory usage
if [ "$use_record" = T ];then
    module load anaconda2
    psrecord_cmd="psrecord --include-children  --interval 1 --log cpu_mem_${outPrefix}.log"
fi


run="matlab -nodisplay -r STBHmftn"  # type waht you mainly want to run here
eval $psrecord_cmd "$run" 1> >(tee ${PBS_JOBID}.${PBS_JOBNAME}.out) 2> >(tee -a ${PBS_JOBID}.${PBS_JOBNAME}.out)
