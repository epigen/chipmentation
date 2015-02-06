#!/usr/env python

import os
import time
import textwrap
import subprocess
import collections
import cPickle as pickle


class DivideAndSlurm(object):
    """
    DivideAndSlurm is a class to handle a map-reduce style submission of jobs to a Slurm cluster.
    Initialize the class with the data to be split, processed in parallel and returned.

    """
    def __init__(self, data, pickleDir="."):
        super(DivideAndSlurm, self).__init__()

        # check data is iterable
        if type(data) == dict:
            data = data.items()
        self._data = data

        self.name = str(time.time())

        self.pickleDir = os.path.abspath(pickleDir)

        self.is_ready = self._is_ready()

    def __repr__(self):
        return "DivideAndSlurm object " + self.name

    def __str__(self):
        return "DivideAndSlurm object " + self.name


    def _chunks(data, n):
        """ Yield successive n-sized chunks from data.
        """
        for i in xrange(0, len(data), n):
            yield data[i:i+n]


    def _slurmHeader(self, jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
        command = """            #!/bin/bash
            #SBATCH --partition={0}
            #SBATCH --ntasks={1}
            #SBATCH --time={2}

            #SBATCH --cpus-per-task={3}
            #SBATCH --mem-per-cpu={4}
            #SBATCH --nodes={5}

            #SBATCH --job-name={6}
            #SBATCH --output={7}

            #SBATCH --mail-type=end
            #SBATCH --mail-user={8}

            # Activate virtual environment
            source /home/arendeiro/venv/bin/activate

            # Start running the job
            hostname
            date

        """.format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output, userMail)

        return command


    def _slurmFooter(self):
        command = """


            # Deactivate virtual environment
            deactivate

            # Job end
            date

        """

        return command


    def _slurmSubmitJob(self, jobFile):
        """
        Submit command to shell.
        """
        command = "sbatch %s" % jobFile
        p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        return p.communicate()


    def split_data(self, fractions):
        """
        Split self._data in fractions and create pickle objects with them.
        """
        chunkify = lambda lst,n: [lst[i::n] for i in xrange(n)]

        groups = chunkify(self._data, len(self._data) / fractions)
        ids = [self.name + "_" + str(i) for i in xrange(len(groups))]
        files = [self.pickleDir + "/" + ID for ID in ids]
        
        # keep track of groups in self
        self.groups = zip(ids, groups, files)

        # serialize groups
        for i in xrange(len(self.groups)):
            pickle.dump(self.groups[i][1],                  # actual group of objects
                open(self.groups[i][2] + ".pickle", 'wb'),  # group pickle file
                protocol=pickle.HIGHEST_PROTOCOL
            )


    def count_distances(self, bam_file, strand_wise=True, duplicates=True, permute=True, fragment_size=1):
        """
        Add task to be performed with data.
        """
        self.jobs = list()
        log = self.name + "_count_distances.log"

        for i in xrange(len(self.groups)):
            jobFile = self.groups[i][2] + "_count_distances.sh"
            input_pickle = self.groups[i][2] + ".pickle"
            output_pickle = self.groups[i][2] + ".output.pickle"

            # assemble command for job
            task = "python count_distances_parallel.py {0} {1} {2} \\".format(input_pickle, output_pickle, bam_file)

            if strand_wise:
                task += "--strand-wise \\"
            if duplicates:
                task += "--duplicates \\"
            if duplicates:
                task += "--permute \\"
            if fragment_size:
                task += "--fragment-size "

            # assemble job file
            job = self._slurmHeader(self.groups[i][0], log) + task + self._slurmFooter()

            # keep track of jobs and their files
            self.jobs.append((job, jobFile))

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))


    def submit(self):
        """
        Submit slurm jobs with each fraction of data.
        """
        jobIDs = list()
        for job, jobFile in self.jobs:
            output, err = _slurmSubmitJob(jobFile)
            jobIDs.append(re.sub("\D", "", output))
        self.submission_time = time.time()
        self.jobIDs = jobIDs


    def _is_ready(self):
        """
        Check if all submitted jobs have been completed.
        """
        # check if all ids are missing from squeue
        for i in xrange(len(self.jobIDs)):
            p = subprocess.Popen("squeue | grep {0}".format(jobIDs[i]), stdout=subprocess.PIPE, shell=True)
            output, err = p.communicate()
            if output.strip() != "":
                return False

        # check if all output pickles are produced
        outputPickles = [self.groups[i][2] + ".output.pickle" for i in xrange(len(self.groups))]
        for i in outputPickles:
            if not os.path.isfile(i):
                return False

        # if both are yes, save output already
        self.output = self._collect_distances()
        return True


    def _collect_distances(self):
        """
        If self.is_ready(), return joined data.
        """
        if self.is_ready:
            # load all pickles into list
            outputs = [pickle.load(open(self.groups[i][2] + ".output.pickle", 'r')) for i in xrange(len(self.groups))]
            # if all are counters, sum them
            if all([type(outputs[i]) == collections.Counter for i in range(len(outputs))]):
                return sum(outputs)
        else:
            return None



        

