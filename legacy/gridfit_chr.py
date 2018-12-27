#!/usr/bin/env python

import datetime
import sys
import os
import time
import subprocess
from argparse import ArgumentParser

devnull = os.open('/dev/null', os.O_WRONLY)
real_stderr = os.dup(2)

def disable_stderr():
    sys.stderr.flush()
    os.dup2(devnull,2)

def enable_stderr():
    sys.stderr.flush()
    os.dup2(real_stderr,2)

def print_job(j):
    print("Completed: {0}".format(j.jobId))
    print("    hasExited          " + str(j.hasExited))
    print("    hasSignal          " + str(j.hasSignal))
    print("    terminatedSignal   " + str(j.terminatedSignal))
    print("    hasCoreDump        " + str(j.hasCoreDump))
    print("    wasAborted         " + str(j.wasAborted))
    print("    exitStatus         " + str(j.exitStatus))
    print("    resourceUsage      " + str(j.resourceUsage))


# from experience: on orchestra, these syncs fail in practice,
# on both long and short chromosomes.  my current guess is that
# lsf on orchestra flushes its job table to disk every hour.
# if this happens and at least one job is still running (so that
# the sync is still waiting), when the sync finally does return,
# the jobs that were flushed to disk fail with
#    e #4908 [ 12362.58]  * call to lsb_openjobinfo returned with error 1:no matching job found mapped to drmaa error 1040:job does not exist in drms queue.
#    e #4908 [ 12362.59]  * fsd_exc_new(1040,lsb_openjobinfo: no matching job found (1),1)
# and throws an invalidjobexception.
# to work around this limitation, use a no-timeout wait and wake
# up every n seconds to recheck the sync status.
# even with a 5 second wait interval, perhaps the invalidjob
# issue could still happen.
def waitall(jobids, session, wait_interval=5):
    print("waiting for %d jobs." % len(jobids))
    nfailed = 0

    still_running = jobids
    while len(still_running) > 0:
        next_running = []
        for j in still_running:
            try:
                # lsf drmaa (the c library) is configured to do pretty
                # verbose tracing and there is no way to disable this from
                # the python bindings.  this is the ugliest way to handle
                # it, but it usually works.
                disable_stderr()
                retval = session.wait(j, drmaa.Session.TIMEOUT_NO_WAIT)
                enable_stderr()
            except drmaa.ExitTimeoutException:
                # Seems like session.wait throws this even for
                # TIMEOUT_NO_WAIT.  maybe there's a good reason
                next_running.append(j)
                continue

            if not retval.hasExited:
                print("Job %d: got retval but job did not exit" % retval.jobId)
                print_job(retval)
                nfailed = nfailed + 1
            else:
                if retval.exitStatus != 0:
                    print_job(retval)
                    nfailed = nfailed + 1
                else:
                    print("Job {0} completed successfully.".format(retval.jobId))

        still_running = next_running
        time.sleep(wait_interval)

    print("All jobs completed.  Checking error statuses..")
    if nfailed > 0:
        print("%d jobs failed. see log above." % nfailed)


def runstep(session, args, k, gridrange, grids_to_run):
    fmtpre = "%sg%d" % (args.outprefix, k)

    print("STEP %d: computing grid" % k)
    rungridC(session, args.queue, args.bindata,
        args.laplace, fmtpre + ".%d",
        gridrange, grids_to_run, args.points_per_grid,
        args.chunksize, args.adj_seed)

    # It is possible that the cluster file system has not yet updated and
    # discover_reruns is based on the expected files produced by each grid.
    # Try to prevent desyncing from launching new jobs by giving the FS a
    # minute to catch up.
    # XXX: only sleep when using DRMAA/cluster
    if session is not None:
        time.sleep(60)

    steps, grids_to_rerun = \
            discover_reruns([ k ], args.ngrids, args.outprefix)
    for i in range(0, args.retries):
        if k in steps:
            print("detected failed grids (retries=%d): %s" \
                % (args.retries, str(grids_to_rerun)))
            print("retry #%d, grid step %d" % (i+1, k))
            rungridC(session, args.queue, args.bindata,
                args.laplace, fmtpre + ".%d",
                gridrange, grids_to_rerun, args.points_per_grid,
                args.chunksize, args.adj_seed)
            steps, grids_to_rerun = \
                    discover_reruns(steps, args.ngrids, args.outprefix)

    if k in steps:
        raise RuntimeError("exhausted retries (%d), giving up" % args.retries)

    print("STEP %d: combining grids" % k)
    return combine(args.combine, args.outprefix, k, args.ngrids)



# uses the C program gridfit-gauss rather than the R script
# gridrange format: length=4, [amin, amax, bmin, bmax]
def rungridC(session, queue, bindata, prog, fmt, gridrange, grids_to_run, ppg, chunksize, adj_seed):
    """Responsible for deciding if jobs should be submitted to a DRM
       (via DRMAA) or run locally."""

    jts = []
    for i in grids_to_run:
        cmd = [ "stdbuf" ] + \
            [ "-o0", "-e0", prog ] + \
            [ bindata, fmt.replace("%d", str(i)) + ".bin" ] +  \
            [ str(x) for x in gridrange ] + \
            [ str(ppg), str(i + adj_seed), str(chunksize) ]

        if session is None:  # running locally
            outputPath = fmt.replace("%d", str(i)) + '.local_log'
            jts.append(' '.join(cmd) + " > " + outputPath)
        else:
            jt = session.createJobTemplate()
            # This disables output buffering to make logs immediately readable
            jt.remoteCommand = cmd[0]
            jt.args = cmd[1:len(cmd)]
            jt.outputPath = ":" + fmt.replace("%d", str(i)) + '.slurm'
            jt.jobName = fmt.replace("%d", str(i)) + '.bin'
            print(jt.args)
            jt.joinFiles = True  # combines stderr and stdout
            queuespec = "-p " + queue
            # submitting to the park queue on O2 requires the park_contrib
            # resource account.
            if queue == 'park':
                queuespec = queuespec + " -A park_contrib"
            jt.nativeSpecification = queuespec + " -t 12:00:00"
            jts.append(jt)
    
    if session is None:  # running locally
        procs = [ subprocess.Popen(jt, shell=True) for jt in jts ]
        for p in procs:
            p.wait()
    else:
        # Doesn't set SLURM_ARRAY_TASK_ID properly
        # I think drmaa's .PARAMETRIC_INDEX might work, but why bother
        #jobids = session.runBulkJobs(jt, 1, ngrids, 1)
        jobids = [ session.runJob(jt) for jt in jts ]
        waitall(jobids, session)
        for jt in jts:
            session.deleteJobTemplate(jt)


def combine(combinescript, outputpre, k, ngrids):
    """Runs an Rscript locally to combine output grids."""
    cmd = [ combinescript, outputpre, str(k), str(ngrids) ]
    output = subprocess.check_output(cmd)
    vals = [ float(line) for line in output.split() ]
    return vals


def combine_and_write(combinescript, outputpre, k, ngrids):
    """Runs an Rscript locally to combine output grids."""
    cmd = [ combinescript, outputpre, str(k), str(ngrids), outputpre + "/fit.rda"]
    output = subprocess.check_output(cmd)
    return True

def discover_reruns(steps, ngrids, outprefix):
    steps_needed = steps[:] # need a copy of steps since we're modifying it
    for stepnum in steps:
        # Each step should have exactly ngrids files named
        # outprefix/g[step].[grid_index].bin
        # A grid needs to be rerun if this binary file is size 0, which
        # will also means this step and all following steps did not finish.
        grids_to_rerun = []
        for i in range(1, ngrids+1):
            fname = "%sg%d.%d.bin" % (outprefix, stepnum, i)
            if not os.path.exists(fname) or os.stat(fname).st_size == 0:
                grids_to_rerun.append(i)

        # The first step with any failed grid is where the resume begins.
        # Steps with all grids intact need not be rerun.
        if len(grids_to_rerun) == 0:
            steps_needed.remove(stepnum)
        else:
            break

    # now steps_needed contains only the steps needed to rerun and
    # grids_to_rerun contains all the failed grids in the last failed step
    return (steps_needed, grids_to_rerun)



def get_gridrange(steps, ngrids, combinescript, outprefix,
                  gridrange=[-7, 2, 2, 4, -7, 2, 2, 6]):
    """determine the grid boundaries for a step. if step is 1, then the
       default range is returned. otherwise, the appropriate grid size
       is calculated from the previous step's results."""
    if len(steps) > 0 and steps[0] > 1:
        gridrange = combine(combinescript, outprefix, steps[0] - 1, ngrids)

    return gridrange



ap = ArgumentParser()
ap.add_argument("--bindata", required=True, metavar="PATH",
    help="Training het sites in binary format from utils/binio.R.")
ap.add_argument("--outprefix", default="", metavar="STRING",
    help="Prepend STRING to all output files.  Can include paths.")
ap.add_argument("--laplace", metavar="PATH",
    default="laplace_cpu",
    help="Path to the binary laplace_cpu approximation program.")
ap.add_argument("--combine", metavar="PATH",
    default="combine_grids.R",
    help="Path to the combine_grids.R script.")
ap.add_argument("--queue", default="park", metavar="STRING",
    help="Submit to SLURM queue STRING.")
ap.add_argument("--ngrids", default=20, metavar="INT", type=int,
    help="Run INT parallel laplace approximators. NOTE: the most recent version samples parameter space randomly rather than separating a rigid grid of points.")
ap.add_argument("--points-per-grid", default=1000, metavar='INT', type=int,
    help='Calculate the likelihood of INT random points within the grid range per grid.')
ap.add_argument("--resume", default=False, action='store_true',
    help="Resume a previous run. Attempts to automatically discover all failed grids, rerun them, and resume the next step if there is one. NOTE: I believe it is always safe to run in resume mode, even when starting a new analysis. One reason to not use resume mode is if you really want to recompute grids that have already completed.")
ap.add_argument("--retries", default=3, metavar="INT", type=int,
    help="When a grid failure is detected, try to rerun that grid INT times before giving up.")
ap.add_argument("--local", default=False, action='store_true',
    help='Run grid search on the local machine; do not submit to a cluster through DRMAA. IMPORTANT: when --local is specified, --ngrids threads will be run in parallel.')
ap.add_argument("--chunksize", default=250, metavar="INT", type=int,
    help='Group hSNPs into non-overlapping windows containing INT hSNPs. Larger chunksizes produce more accurate parameter estimates; however, runtime increases by O(N^3).')
ap.add_argument("--adj-seed", default=0, metavar="INT", type=int,
    help="Change the random seed used for selecting random points in (a,b,c,d) parameter space. By default, the seed is equal to the grid ID.  I.e., if parallelizing over 20 grids, grid N would be given seed N + INT (default: INT=0). This option is intended for testing the stability of parameter estimation and is not of practical use.")
args = ap.parse_args()

args.outprefix = os.path.realpath(args.outprefix) + '/'

s = None
if args.local == False:
    import drmaa
    s = drmaa.Session()
    s.initialize()
    print('Supported contact strings: %s' % s.contact)
    print('Supported DRM systems: %s' % s.drmsInfo)
    print('Supported DRMAA implementations: %s' % s.drmaaImplementation)

steps = [ 1, 2, 3, 4 ]
grids_to_run = range(1, args.ngrids+1)
if args.resume:
    print("Resuming")
    steps, grids_to_run = \
        discover_reruns(steps, args.ngrids, args.outprefix)

    if len(steps) > 0:
        print("Failed grids from step %d are %s" % \
            (steps[0], str(grids_to_run)))

gridrange = get_gridrange(steps, args.ngrids, args.combine, args.outprefix)
print("Running steps " + str(steps))
print("Initial grid ranges " + str(gridrange))

# use an objective stopping criteria one day
# for now, empirically determined that parameters tend to be very
# converged by g3, and g4 is overkill.
print('Starting run at ' + str(datetime.datetime.now()))
start_time = time.time()
for stepnum in steps:
    newgrid = runstep(s, args, stepnum, gridrange, grids_to_run)
    print(newgrid)
    gridrange = newgrid[0:8]

if args.local == False:
    s.exit()

#def combine_and_write(combinescript, outputpre, k, ngrids):
print("saving best fit from grid " + str(steps[len(steps)-1]))
combine_and_write(args.combine, args.outprefix, steps[len(steps)-1], args.ngrids)

print('Run completed ' + str(datetime.datetime.now()))
end_time = time.time()
print("%0.2f total seconds elapsed" % (end_time - start_time))
