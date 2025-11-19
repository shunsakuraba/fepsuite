# Overview

Files in this directory are scripts for supercomputer job submission systems. There are no 'standard' method of the job submission - even if two supercomputer center use the same job system, they often have customization specific to one supercomputer center and they almost never share the same MPI calling commands.
The files here are my solution to this mess - I try to absorb the difference among centers as possible, with a bunch of scripting files.

# What you need to do as a user

If you are not a user of a supercomputer system, the best method you should choose is run with "`JOBSYSTEM=none`", where everything is assumed to finish within your machine. We also assume the toolchains and environments follow current standards - OpenMPI is installed and the node has one or more GPUs.

If you want to use a supercomputer, and if you can find corresponding scripts, specify that to `JOBSYSTEM`. You may probably need to fill in "`GROUP`" and/or "`SUBSYSTEM`" variables to use the script (read the error message if you need).

If you want to use a supercomputer, and if you can't find the script, then you need to proceed to the next section.

# What you need to do as a developer

You need to implement 8 shell functions. Open template.zsh and you can find more detailed information.

`job_prelude()`: Change directory to submitted one and load modules. You also need to set PYTHON3 to be python3

`job_prelude_after_gmx()`: Probably you should do nothing here. If you need post-hook after source GMXRC write something here.

`job_init_queue_stat()`: Probably you should do nothing here.

`job_set_preferred_resource()`: This function typically just sets `CPN` (preferred core per node), `GPN` (GPU per node), `HW_CPN` (REAL core per node, typically =`CPN`), and `HW_GPN` (REAL GPU per node). `GPN` value may be greater than `HW_GPN`; e.g., even if your machine may have only 2 GPUs, you may want to run 4 or 8 GPU-using processes for a better *throughput*. You can also change the resource limitation based on requested resources - this is useful when e.g. the supercomputer have a dedicated smaller resource queue (e.g. "shared node").

`job_submit()`: This function submits jobs based on requested resources. You need to consider four things: (1) export necessary environment variables to the job file, (2) add dependency condition to job (or wait until the job finishes), (3) write proper resource requirement to tell job system, and (4) return job ID which will be used in the future job submission as "dependency" info.

`job_mpirun()`: Run MPI with the number of ranks with the first argument (`$1`), number of threads `$OMP_NUM_THREADS`, and with commands specified in `$@` (you need to remove `$1` by `shift` command first). If your MPI caller can control thread pinning, use `$PPN` as processes per node.

`job_singlerun()`: Run program specified in arguments; this is typically just running `$@`, but some supercomputers allow you to only run programs via specific commands (`srun -n 1` , `aprun -n 1`, ...)

`job_get_mode()`: This script is called in two ways - one when you submit the job, and the other when you run the batch job. You must return whether it's inside job or not (e.g. by checking environment variables exists)

# Concepts

In these submission scripts, resources are abstracted so as to follow the terminology.

<dl>
  <dt>Cores</dt>
  <dd>Typically refers to the number of physical CPU cores. Because using the same core by two different CPU tasks significantly slows down the calculation (and hyper-threading is rarely useful in this field), number of cores in the node (also see "Nodes" below) is used as the upper limit of the number of threads.</dd>
  <dt>Processes</dt>
  <dd>The number of MPI processes. Processes are decided by two factors: number of processes to run each simulation (e.g. `*_PARA=` in abfe and `PARA=` in feprest), and number of simulations simultaneously run by the replica exchange. </dd>
  <dt>Threads</dt>
  <dd>The number of CPU threads per each process. In pure CPU run it is recommended to be 1. In GPU run it is recommended to be around 2-8.</dd>
  <dt>Nodes</dt>
  <dd>The "node" is a set of resource that is *not* splittable in your job system. Cores less than one node are round up so that it is the multiple of `CPN`. Based on CPN/GPN and required MPI processes, the script adjusts the *number of nodes* upon the job submission. Calculated number of nodes are represented as `JOB_NODES`. Because some supercomputer systems only allow us to set the number of processes but /not number of nodes/, the script can use `JOB_PROC` variable, which is number of process calculated from number of nodes upon submission.</dd>
</dl>

# About bundled scripts

You may want to use the script below as an alternative template.

<dl>
  <dt>nagoya-cx.zsh</dt>
  <dd>Nagoya University, Fujitsu's job manager (pjsub), OpenMPI</dd>
  <dt>kudpc-2023.zsh</dt>
  <dd>Kyoto University, slurm (but their resource request part is customized), Intel MPI (MPICH-like)</dd>
  <dt>jaea-sgi8600.zsh</dt>
  <dd>Japan Atomic Energy Agency, PBSPro, Intel MPI + SGI mpt</dd>
  <dt>tokyo-wisteria.zsh</dt>
  <dd>The University of Tokyo, Wisteria, Fujitsu pjsub, Intel MPI</dd>
  <dt>kyushu-genkai.zsh</dt>
  <dd>Kyushu University, Genkai, Fujitsu pjsub, Intel MPI</dd>
  <dt>tsukuba-pegasus.zsh</dt>
  <dd>Tsukuba Univresity, Pegasus, NQSV, OpenMPI</dd>
</dl>

# Tips for zsh scripting

I use zsh for writing scripts because it allows me to calculate (not so) complex maths inside scripts. But some of you may be experienced user in bash and not in zsh. For them, here are some tips.

Zsh's array index starts from 1. `${ARR[1]}` is the first element. The number of elements is `${#ARR}`.

Zsh's command invocation does not automatically split the command line by " ". For example:
````sh
foo="X Y"
tr $foo <<< "Input XYZ" # tr: missing operand after ‘X Y’
````
will be *an error*. This is because `"X Y"` is *directly* sent to `tr` as the second argument, without splitting by `" "`. Instead you need to use arrays:
````sh
foo=(X Y)
tr $foo <<< "Input XYZ"
````

Expanding variable from variable names can be done with ${(P)x}.

````sh
foo=var
var=1
echo ${(P)foo} # "1"
````

To join an array with a character:
````sh
foo=(X Y Z)
echo ${(j:,:)foo} # "X,Y,Z"
echo ${(j|:|)foo} # "X:Y:Z"
````

To split string into array:
````sh
foo="X,Y,Z"
ARR=(${(s:,:)foo}) # split by ','
echo ${ARR[2]}     # "Y"
````


