#!/usr/env python

"""
Predicting nucleosome associated-regions from ChIPmentation data.

In brief, this does:
    * Pairwise distance count between reads from several samples (from several
    techniques - ChIP, ChIPmentation, DNase-seq, ATAC-seq) in several genomic
    locations (H3K4me3/H3K27me3 peaks, DHS regions, whole genome, etc...) in strand-
    -specific way and also for permuted reads.
    * Get pattern from the read distances distribution (correlogram).
    * Decompose pattern into signals with various frequencies (using FFT) and retrieve
    most abundant (10bp).
    * Calculate correlations between pattern and ChIPmentation read coverage/permuted reads along the genome.
    * Extract local maxima from correlation into a binary signal.
    * Feed to a HMM modeling a nucleosome, output predicted nucleosome-associated regions.
    * Measure several features of predicted regions.
    * Plot and calculate p-values for each feature between real data and permuted.
    * Train classifier with features from real and permuted data.
    * Classify the rest of data N times independently and calculate FDR for each nucleosome-associated region.
    * Export nucleosome-associated regions with FDR<0.5.

    * ... plus several plots along the way.

TODO:
    * Plot CM, DNase, MNase in DARNS:
        * From the middle of the DARN (mid-peak)
        * From the 5' and 3' end of nucleosome Dyads (get external data)
        * DARNS frequency around TSSs (models and CAGE) and TTSs
        * DARNS frequency around CpGs islands

    * Try:
        * Do the same on concatenated histone ChIPmentation data.
        * Do the same using IGG data as background.
"""

from collections import Counter, OrderedDict
from divideAndSlurm import DivideAndSlurm, Task
from matplotlib import pyplot as plt
from pybedtools import BedTool
import cPickle as pickle
import HTSeq
import itertools
import numpy as np
import pandas as pd
import os
import random
import re
import seaborn as sns
import string
import textwrap
import time
import yahmm
from sklearn.metrics import classification_report, roc_curve, roc_auc_score, precision_recall_curve, average_precision_score

sns.set_style("whitegrid")
np.set_printoptions(linewidth=200)


class CountDistances(Task):
    """
    Task to perform counting of distances between reads under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(CountDistances, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "count_distances_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Bam file argument is missing")
        self.bam_file = self.args[0]
        # additional arguments
        if "strand_wise" in kwargs.keys():
            self.strand_wise = kwargs["strand_wise"]
        else:
            self.strand_wise = True
        if "duplicates" in kwargs.keys():
            self.duplicates = kwargs["duplicates"]
        else:
            self.duplicates = True
        if "permute" in kwargs.keys():
            self.permute = kwargs["permute"]
        else:
            self.permute = False
        if "fragment_size" in kwargs.keys():
            self.fragment_size = kwargs["fragment_size"]
        else:
            self.fragment_size = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_count_distances.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python count_distances_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.permute:
                task += "--permute "
            task += "--fragment-size {0}".format(self.fragment_size)

            task += """

                # Deactivate virtual environment
                deactivate
                    """

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are counters, and their elements are counters, sum them
                if all([type(outputs[i]) == Counter for i in range(len(outputs))]):
                    output = reduce(lambda x, y: x + y, outputs)  # reduce
                    if type(output) == Counter:
                        self.output = output    # store output in object
                        self._rm_temps()  # delete tmp files
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class Coverage(Task):
    """
    Task to get read coverage under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(Coverage, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "coverage_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Bam file argument is missing")
        self.bam_file = self.args[0]
        # additional arguments
        if "strand_wise" in kwargs.keys():
            self.strand_wise = kwargs["strand_wise"]
        else:
            self.strand_wise = False
        if "duplicates" in kwargs.keys():
            self.duplicates = kwargs["duplicates"]
        else:
            self.duplicates = True
        if "orientation" in kwargs.keys():
            self.orientation = kwargs["orientation"]
        else:
            self.orientation = False
        if "permute" in kwargs.keys():
            self.permute = kwargs["permute"]
        else:
            self.permute = False
        if "fragment_size" in kwargs.keys():
            self.fragment_size = kwargs["fragment_size"]
        else:
            self.fragment_size = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_coverage.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python coverage_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.orientation:
                task += "--orientation "
            if self.permute:
                task += "--permute "
            task += "--fragment-size {0}".format(self.fragment_size)

            task += """

                # Deactivate virtual environment
                deactivate
                    """

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are counters, and their elements are counters, sum them
                if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                    output = reduce(lambda x, y: dict(x, **y), outputs)
                    if type(output) == dict:
                        self.output = output  # store output in object
                        self._rm_temps()  # delete tmp files
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class Correlation(Task):
    """
    Task to get read coverage under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(Correlation, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "correlation_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Pattern argument is missing")
        self.pattern = self.args[0]
        # additional arguments
        if "step" in kwargs.keys():
            self.step = kwargs["step"]
        else:
            self.step = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Pickle pattern
        self.patternPickle = os.path.join(self.slurm.tmpDir, self.name + "_pattern.pickle")
        pickle.dump(self.pattern, open(self.patternPickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_correlation.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath to python script!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python correlation_parallel.py {0} {1} {2} --step {3}

                # Deactivate virtual environment
                deactivate
                """.format(inputPickle, outputPickle, self.patternPickle, self.step)

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are dicts, and their elements are dicts, sum them
                if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                    output = reduce(lambda x, y: dict(x, **y), outputs)
                    if type(output) == dict:
                        self.output = output  # store output in object
                        self._rm_temps()  # delete tmp files
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class CorrelatePatternBam(Task):
    """
    Task to get Correlation peaks from reads and a given pattern.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(CorrelatePatternBam, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "CorrelatePatternBam_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 2:
            raise TypeError("You must specify two arguments: Pattern, Bam")
        self.pattern = self.args[0]
        self.bam_file = self.args[1]
        # additional arguments
        if "strand_wise" in kwargs.keys():
            self.strand_wise = kwargs["strand_wise"]
        else:
            self.strand_wise = True
        if "duplicates" in kwargs.keys():
            self.duplicates = kwargs["duplicates"]
        else:
            self.duplicates = True
        if "orientation" in kwargs.keys():
            self.orientation = kwargs["orientation"]
        else:
            self.orientation = False
        if "permute" in kwargs.keys():
            self.permute = kwargs["permute"]
        else:
            self.permute = False
        if "fragment_size" in kwargs.keys():
            self.fragment_size = kwargs["fragment_size"]
        else:
            self.fragment_size = 1
        if "step" in kwargs.keys():
            self.step = kwargs["step"]
        else:
            self.step = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Pickle pattern
        self.patternPickle = os.path.join(self.slurm.tmpDir, self.name + "_pattern.pickle")
        pickle.dump(self.pattern, open(self.patternPickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_correlatePatternBam.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath to python script!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python correlatePatternBam_parallel.py {0} {1} {2} {3} """.format(inputPickle, outputPickle, self.patternPickle, self.bam_file)

            task += "--step {0} ".format(self.step)
            task += "--fragment-size {0} ".format(self.fragment_size)
            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.orientation:
                task += "--orientation "
            if self.permute:
                task += "--permute "

            task += """

                # Deactivate virtual environment
                deactivate
                    """

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are dicts, and their elements are dicts, sum them
                if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                    output = reduce(lambda x, y: dict(x, **y), outputs)
                    if type(output) == dict:
                        # Split name into chrom start end
                        output = {(key.split("_")[0], int(key.split("_")[1]), int(key.split("_")[2])): value for key, value in output.items()}
                        # Sort, make ordered dict, store output in object
                        self.output = OrderedDict(sorted(output.items(), key=lambda x: (x[0][0], x[0][1])))
                        # Delete tmp files
                        self._rm_temps()
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class WinterHMM(object):
    """
    Class to model the Hidden Markov model described in the Winter 2013 paper.
    """
    def __init__(self):
        super(WinterHMM, self).__init__()

        self.mapping = {v: k for k, v in enumerate(itertools.product([1, 0, -1], repeat=2))}
        # self.mapping = {
        #     (-1, -1):   8,
        #     (-1, 0):    7,
        #     (-1, 1):    6,
        #     (0, -1):    5,
        #     (0, 0):     4,
        #     (0, 1):     3,
        #     (1, -1):    2,
        #     (1, 0):     1,
        #     (1, 1):     0
        # }

        # rewrite like this:
        # matrix = array([
        #     [-inf,-inf, -0.69314718,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf, -0.69314718,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.10536052,-inf,-inf,-inf,-inf, -2.30258509, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718, -0.69314718,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -0.7985077 ,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.7985077 , -2.30258509, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf, -inf,-inf],
        #     [-0.69314718,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -2.30258509,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -2.30258509,-inf,-inf, -0.22314355, -inf,-inf],
        #     [-inf,-inf, -2.30258509,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -2.30258509,-inf,-inf, -0.22314355, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf]
        # ])

        # distributions = [NormalDistribution(1, .5), NormalDistribution(5, 2)]
        # starts = [ 1., 0. ]
        # ends = [ .1., .1 ]
        # state_names= [ "A", "B" ]

        # model = Model.from_matrix( matrix, distributions, starts, ends,
        #         state_names, name="test_model" )

        background = yahmm.DiscreteDistribution({
            0: 0.001, 1: 0.001, 2: 0.001,
            3: 0.001, 4: 0.99, 5: 0.002,
            6: 0.001, 7: 0.002, 8: 0.001
        })
        pos = yahmm.DiscreteDistribution({
            0: 0.33, 1: 0.33, 2: 0.33,
            3: 0.001, 4: 0.001, 5: 0.001,
            6: 0.001, 7: 0.001, 8: 0.001
        })
        neg = yahmm.DiscreteDistribution({
            0: 0.33, 1: 0.001, 2: 0.001,
            3: 0.33, 4: 0.001, 5: 0.001,
            6: 0.33, 7: 0.001, 8: 0.001
        })
        zero = yahmm.DiscreteDistribution({
            0: 0.001, 1: 0.2, 2: 0.001,
            3: 0.2, 4: 0.2, 5: 0.2,
            6: 0.001, 7: 0.2, 8: 0.001
        })

        # Create states
        b = yahmm.State(background, name="B")
        minus = yahmm.State(neg, name="-")
        in1 = yahmm.State(zero, name="in1")
        in2 = yahmm.State(zero, name="in2")
        in3 = yahmm.State(zero, name="in3")
        plus = yahmm.State(pos, name="+")
        in4 = yahmm.State(zero, name="in4")
        in5 = yahmm.State(zero, name="in5")
        in6 = yahmm.State(zero, name="in6")
        in7 = yahmm.State(zero, name="in7")
        in8 = yahmm.State(zero, name="in8")
        in9 = yahmm.State(zero, name="in9")
        in10 = yahmm.State(zero, name="in10")
        in11 = yahmm.State(zero, name="in11")

        self.model = yahmm.Model(name="winter-2013")

        for state in [b, minus, in1, in2, in3, plus, in4, in5, in6, in7, in8, in9, in10, in11]:
            self.model.add_state(state)

        self.model.add_transition(self.model.start, b, 0.8)       # start in background
        self.model.add_transition(self.model.start, minus, 0.1)   # start in minus
        self.model.add_transition(self.model.start, plus, 0.1)    # start in plus

        self.model.add_transition(b, b, 0.8)                      # stay in background
        self.model.add_transition(b, minus, 0.1)                 # enter in minus
        self.model.add_transition(b, plus, 0.1)                  # enter in plus

        self.model.add_transition(minus, b, 0.1)                  # can exit from plus
        self.model.add_transition(minus, in1, 0.45)
        self.model.add_transition(minus, plus, 0.45)

        self.model.add_transition(in1, in2, 0.5)                  # states 1-3 can go directly to plus or follow a string 1 to 4
        self.model.add_transition(in1, plus, 0.5)
        self.model.add_transition(in2, in3, 0.5)
        self.model.add_transition(in2, plus, 0.5)
        self.model.add_transition(in3, in4, 0.5)
        self.model.add_transition(in3, plus, 0.5)

        self.model.add_transition(plus, b, 0.1)                   # can exit from plus
        self.model.add_transition(plus, in4, 0.9)
        self.model.add_transition(in4, in5, 1)                    # string of 4 to 8 at least
        self.model.add_transition(in5, in6, 1)
        self.model.add_transition(in6, in7, 1)
        self.model.add_transition(in7, in8, 1)
        self.model.add_transition(in8, in9, 0.5)                  # states 8-11 can go directly to minus
        self.model.add_transition(in8, minus, 0.5)
        self.model.add_transition(in9, in10, 0.5)
        self.model.add_transition(in9, minus, 0.5)
        self.model.add_transition(in10, in11, 0.5)
        self.model.add_transition(in10, minus, 0.5)
        self.model.add_transition(in11, minus, 1)

        self.model.bake()

    def draw(self):
        """
        Draws the Markov chain of the model.
        """
        return self.model.draw(node_size=400, labels={state.name: str(state.name) for state in self.model.states}, font_size=20)

    def train(self, sequences):
        """
        Train model with given data.
        """
        self.model.train(sequences)

    def predict(self, observations):
        """
        Predict hidden states from observations.
        """
        logp, chain = self.model.maximum_a_posteriori(observations)
        return [state.name for num, state in chain][1:-1]  # strip begin and end states

    def retrieveProbabilities(self, observations):
        """
        Retrieve the posterior probabilities of being in a "nucleosome" state for a sequence of observations.
        """
        trans, ems = self.model.forward_backward(observations)
        ems = np.exp(ems)
        self.probs = ems / np.sum(ems, axis=1)[:, np.newaxis]  # probs of all states
        background_state = 0
        prob_nucleosome = 1 - self.probs[:, background_state]  # prob of all states but background
        return prob_nucleosome

    def retrieveBackgroundProbabilities(self, observations):
        """
        Retrieve the posterior probabilities of being in a "nucleosome" state for a sequence of observations.
        """
        trans, ems = self.model.forward_backward(observations)
        ems = np.exp(ems)
        self.probs = ems / np.sum(ems, axis=1)[:, np.newaxis]  # probs of all states
        background_state = 0
        prob_background = self.probs[:, background_state]  # prob of all states but background
        return prob_background


class MeasureDarnsFeatures(Task):
    """
    Task to measure features of predicted DARNS.
    Takes as arguments an iterable with indexes of darns in a specific chromosome and type of DARNS (real or permuted)
    and returns dataframe with their measurements.
    """
    def __init__(self, data, fractions, *args):
        super(MeasureDarnsFeatures, self).__init__(data, fractions, *args)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "MeasureDarnsFeatures_{0}".format(now)

        # Parse
        # required arguments
        if len(self.args) != 2:
            raise TypeError("Invalid number of arguments passed to creat MeasureDarnsFeatures object.")
        self.regionType = self.args[0]
        self.chrom = self.args[1]

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_MeasureDarnsFeatures.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath!
            task = """

                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python measure_darns_features_parallel.py {0} {1} {2} {3}

                # Deactivate virtual environment
                deactivate
                """.format(inputPickle, outputPickle, self.regionType, self.chrom)

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are counters, and their elements are counters, sum them
                if all([isinstance(outputs[i], type(pd.Series())) for i in range(len(outputs))]):
                    output = pd.DataFrame([isinstance(outputs[i], type(pd.Series())) for i in range(len(outputs))])  # join in dataframe

                    if isinstance(output, type(pd.DataFrame())):
                        self.output = output    # store output in object
                        self._rm_temps()  # delete tmp files
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


def measureDarnsFeatures(regions, hmm, chrom, regionType, features, coverage, darnsNumber):
    """
    """
    if darnsNumber % 10000 == 0:
        print darnsNumber

    start, end, center = regions[chrom][darnsNumber]
    name = "_".join([chrom, str(start), str(end), str(center)])

    series = pd.Series(index=["type"] + ["name"] + features)
    series["type"] = regionType
    series["name"] = name

    # measure length
    series["length"] = end - start

    # measure distance to neighbours
    if darnsNumber != 0:
        series["space_upstream"] = abs(start - regions[chrom][darnsNumber - 1][1])  # current start minus previous end
    else:
        series["space_upstream"] = None
    if darnsNumber != len(regions[chrom]) - 1:
        series["space_downstream"] = abs(regions[chrom][darnsNumber + 1][0] - end)  # current end minus next start
    else:
        series["space_downstream"] = None

    # get posterior prob
    # sequence = genome_binary[chrom][start: end]
    # series["post_prob"] = model.retrieveProbabilities(sequence)

    # get read count
    series["read_count"] = coverage[regionType][name].sum()
    # get density
    # I defined it as (sum / length). It is not clear this is the same as in Winter2013
    series["read_density"] = coverage[regionType][name].sum() / (end - start)

    # n. of positive peaks in nucleosome model
    hmm = eval("hmmOutput") if regionType == "DARNS" else eval("hmmOutputP")
    series["n_pospeaks"] = hmm[chrom][start: end].count("+")

    # n. of negative peaks in nucleosome model
    series["n_negpeaks"] = hmm[chrom][start: end].count("-")

    # append
    return series


def makeGenomeWindows(windowWidth, genome, step=None):
    """
    Generate windows genome-wide for a given genome with width=windowWidth and
    return dictionary of HTSeq.GenomicInterval objects.

    windowWidth=int.
    genome=str.
    """
    if step is None:
        step = windowWidth
    w = BedTool.window_maker(BedTool(), genome=genome, w=windowWidth, s=step)
    windows = dict()
    for interval in w:
        feature = HTSeq.GenomicInterval(
            interval.chrom,
            interval.start,
            interval.end
        )
        name = string.join(interval.fields, sep="_")
        windows[name] = feature

    return windows


def makeBedWindows(windowWidth, bedtool, step=None):
    """
    Generate windows with width=windowWidth given a pybedtools.BedTool object and
    return dictionary of HTSeq.GenomicInterval objects.

    windowWidth=int.
    bedtool=pybedtools.BedTool object.
    """
    if step is None:
        step = windowWidth
    w = BedTool.window_maker(BedTool(), b=bedtool, w=windowWidth, s=step)
    windows = dict()
    for interval in w:
        feature = HTSeq.GenomicInterval(
            interval.chrom,
            interval.start,
            interval.end
        )
        name = string.join(interval.fields, sep="_")
        windows[name] = feature

    return windows


def bedToolsInterval2GenomicInterval(bedtool, strand=True, name=True):
    """
    Given a pybedtools.BedTool object, return dictionary of HTSeq.GenomicInterval objects.

    bedtool=pybedtools.BedTool object.
    """
    intervals = OrderedDict()
    if strand:
        for iv in bedtool:
            if name:
                intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
            else:
                intervals[string.join(iv.fields[:3], sep="_")] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
    else:
        for iv in bedtool:
            if name:
                intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)
            else:
                intervals[string.join(iv.fields[:3], sep="_")] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)

    return intervals


def DARNS2GenomicInterval(DARNS):
    """
    """
    intervals = OrderedDict()
    for chrom, values in DARNS.items():
        for start, end, center in values:
            name = "_".join([chrom, str(start), str(end), str(center)])
            intervals[name] = HTSeq.GenomicInterval(chrom, start, end)

    return intervals


def extractPattern(dists, distRange, filePrefix):
    """
    Performs the fast Fourier transform (fft) on a dict of x:y values in several frequencies, and
    returns the signal of the most abundant (highest signal amplitude) frequency through inverse fft.
    Produces plots of the process.

    dists=dict - distances:counts, both int.
    distRange=int - indexes of subset of dists to extract pattern from.
    filePrefix=str - prefix of files to save plots to (should include paths).
    """
    plt.close()
    plt.figure(0)
    plt.subplot(211)
    x = dists.keys()
    y = [dists[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1))  # fit linear regression
    m, b = p1.coeffs
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )

    # restrict to signal
    x = distRange
    y = [dists[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1))  # fit linear regression
    m, b = p1.coeffs
    plt.subplot(212)
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )
    plt.savefig(filePrefix + ".read_distances.pdf", bbox_inches='tight')
    plt.close()

    # measure distance to regression in that point
    distsReg = [y[i] - (m * i + b) for i in range(len(x))]

    # subtract value of minimum to all points
    distsReg -= min(distsReg)

    # fourier transform data
    time = x
    signal = distsReg

    # get frequencies from decomposed fft
    W = np.fft.fftfreq(signal.size, d=time[1] - time[0])
    f_signal = np.fft.fft(signal)

    # plot all frequencies, ask what is the amplitude of the signals
    freqs = dict()
    plt.figure(1)
    for i in np.arange(0, len(x) / 10., 0.01):
        if i == 0:
            continue
        else:
            cut_f_signal = f_signal.copy()
            cut_f_signal[((W < i) | (W > i))] = 0
            cut_signal = np.fft.ifft(cut_f_signal)
            plt.plot(time, cut_signal)
            if 0.05 <= i <= 0.4:
                freqs[i] = np.abs(cut_f_signal).max()
    plt.savefig(filePrefix + ".fft_frequencies.pdf", bbox_inches='tight')
    plt.close()

    # get frequency of signal with highest amplitude
    top = max(freqs, key=freqs.get)

    # signal is now in Hz
    cut_f_signal = f_signal.copy()

    # select frequency of top/10bp
    cut_f_signal[((W < top) | (W > top))] = 0

    # inverse fourier to get filtered frequency
    cut_signal = np.fft.ifft(cut_f_signal)

    plt.figure(2)
    plt.subplot(221)
    plt.plot(time, signal, '-')
    plt.subplot(222)
    plt.plot(W, abs(f_signal), 'o')
    plt.subplot(223)
    plt.plot(W, abs(cut_f_signal), 'o')
    plt.subplot(224)
    plt.plot(time, cut_signal, '-')
    plt.savefig(filePrefix + ".fft_filter-{0}bp_ifft.pdf".format(int(1 / top)), bbox_inches='tight')
    plt.close()

    return cut_signal.real


def concatenateBinary(binaryDict, patternLength, debug=False):
    """
    Concatenates string of binary signals from several consecutive windows.
    Returns dict of chr:binary.

    binaryDict=dict - name: 1D numpy.array of ~binary (1,0,-1) signals.
    patternLength=int - length of pattern used to create binary correlation peaks.
    """
    # TODO: make sure windows are sorted alphanumerically
    items = binaryDict.items()
    binaries = OrderedDict()
    prev_chrom = ""
    prev_end = 1

    for i in xrange(len(items)):
        window, sequence = items[i]
        chrom, start, end = (str(window[0]), int(window[1]), int(window[2]))
        if debug:
            name = "_".join([str(j) for j in window])
            print name, (chrom, start, end)
        if not chrom == prev_chrom:  # new chromossome
            binaries[chrom] = sequence
        else:  # same chromossome as before
            if window == items[-1][0]:  # if last window just append remaining
                if debug:
                    print("Last of all")
                binaries[chrom] = np.hstack((binaries[chrom], sequence))
            elif items[i + 1][0][0] != chrom:  # if last window in chromossome, just append remaining
                if debug:
                    print("Last of chromossome")
                binaries[chrom] = np.hstack((binaries[chrom], sequence))
            else:  # not last window
                if prev_end - patternLength == start:  # windows are continuous
                    if len(sequence) == (end - start) - patternLength + 1:
                        binaries[chrom] = np.hstack((binaries[chrom], sequence))
                    elif len(sequence) > (end - start) - patternLength + 1:
                        if debug:
                            print(name, len(sequence), (end - start) - patternLength + 1)
                        raise ValueError("Sequence is bigger than its coordinates.")
                    elif len(sequence) < (end - start) - patternLength + 1:
                        if debug:
                            print(name, len(sequence), (end - start) - patternLength + 1)
                        raise ValueError("Sequence is shorter than its coordinates.")
                else:
                    if debug:
                        print(name, prev_end, start)
                    raise ValueError("Windows are not continuous or some are missing.")

        prev_chrom = chrom
        prev_end = end
    return binaries


def splitAndSortDict(binary):
    """
    Splits dictionary keys in chrom, start, end; sorts dict based on that and returns collections.OrderedDict object.

    binary=dict - keys must be strings that when separated give a chrom (string), start and end (both ints).
    """
    b = {(key.split("_")[0], int(key.split("_")[1]), int(key.split("_")[2])): value for key, value in binary.items()}
    return OrderedDict(sorted(b.items(), key=lambda x: (x[0][0], x[0][1])))


def getDARNS(sequence, debug=False):
    """
    Returns list of tuples with start and end positions (0-based) of DARNS.
    DARNS are stretches of sequences which belong to the "nucleosome" states (see WinterHMM object definition).

    sequence=str.
    """
    outside = ["B"]
    DARNS = list()
    DARN = False
    l = len(sequence)

    # Detect DARNS start and end
    for i in xrange(l):
        if not i == l - 1:  # not last element
            if debug: print(i, "Not last element")
            # Entering DARN
            if not DARN and (sequence[i] not in outside) and i != len(sequence):
                if debug: print(i, "Entered DARN")
                DARN = True
                DARN_start = i
                # Handle 1-element DARNS
                if sequence[i + 1] in outside:
                    if debug: print(i, "1bp DARN")
                    #DARNS.append((DARN_start, i + 1))
                    DARN = False
            # DARN ends at next element
            elif DARN and (sequence[i] not in outside) and (sequence[i + 1] in outside):
                if debug: print(i, "Exited DARN")
                DARN = False
                DARN_end = i
                DARNS.append((DARN_start, DARN_end))
        elif i == l - 1 and DARN:  # last element
            if debug: print(i, "Last element")
            if DARN:  # finish DARN at end
                if debug: print(i, "Finishing DARN in last element")
                DARNS.append((DARN_start, i))

    # Find DARNS middle point
    DARNS_middle = list()
    for start, end in DARNS:
        seq = sequence[start: end + 1]
        middle = int(((end + 1) - start) / 2)
        negative_peaks = [i for i in xrange(len(seq)) if seq[i] == '-']

        if len(negative_peaks) is not 0:
            # find middle most negative peak
            peak = start + min(negative_peaks, key=lambda x: abs(x - middle))
        else:
            peak = start + middle

        DARNS_middle.append((start, end, peak))

    return DARNS_middle


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True, strand_specific=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    stranded - boolean.
    duplicates - boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']
    cov = OrderedDict()
    n = len(intervals)
    i = 0
    for name, feature in intervals.iteritems():
        if i % 1000 == 0:
            print(n - i)
        # Initialize empty array for this feature
        if not strand_specific:
            profile = np.zeros(feature.length, dtype=np.float64)
        else:
            profile = np.zeros((2, feature.length), dtype=np.float64)

        # Check if feature is in bam index
        if feature.chrom not in chroms or feature.chrom == "chrM":
            i += 1
            continue

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            # check it's aligned
            if not aln.aligned:
                continue

            aln.iv.length = fragmentsize  # adjust to size

            # get position in relative to window
            if orientation:
                if feature.strand == "+" or feature.strand == ".":
                    start_in_window = aln.iv.start - feature.start - 1
                    end_in_window = aln.iv.end - feature.start - 1
                else:
                    start_in_window = feature.length - abs(feature.start - aln.iv.end) - 1
                    end_in_window = feature.length - abs(feature.start - aln.iv.start) - 1
            else:
                start_in_window = aln.iv.start - feature.start - 1
                end_in_window = aln.iv.end - feature.start - 1

            # check fragment is within window; this is because of fragmentsize adjustment
            if start_in_window <= 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                profile[start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window: end_in_window] += 1
                else:
                    profile[1][start_in_window: end_in_window] += 1

        # append feature profile to dict
        cov[name] = profile
        i += 1
    return cov


def measureDARNS(DARNS):
    """
    Computes distributions of DARNS attributes: width, interdistance, distance between midpoints.

    tuples=list - list of tuples: (start,end).
    """
    start, end = range(2)
    distances = Counter()
    midDistances = Counter()

    try:
        for chrom, darns in DARNS.items():
            for d1, d2 in itertools.combinations(darns, 2):
                # distance end-to-end
                if d1[end] <= d2[start]:
                    distances[abs(d2[start] - d1[end])] += 1
                else:
                    distances[d1[start] - d2[end]] += 1

                # distance midpoint-midpoint
                midDistances[abs(((d1[end] + d1[start]) / 2) - ((d2[end] + d2[start]) / 2))] += 1
        return (distances, midDistances)
    except KeyboardInterrupt:
        return (distances, midDistances)


def exportWigFile(intervals, profiles, offset, filename, trackname):
    """
    Exports a wig file track with scores contained in profile, .

    intervals=iterable with HTSeq.GenomicInterval objects.
    profiles=iterable with scores to write.
    offset=int.
    filename=str.
    trackname=str.
    """
    with open(filename, 'w') as handle:
        track = 'track type=wiggle_0 name="{0}" description="{0}" visibility=dense autoScale=off\n'.format(trackname)
        handle.write(track)
        for i in xrange(len(profiles)):
            header = "fixedStep  chrom={0}  start={1}  step=1\n".format(intervals[i].chrom, intervals[i].start + offset)
            handle.write(header)
            for j in xrange(len(profiles[i])):
                handle.write(str(abs(profiles[i][j])) + "\n")


def exportBedFile(intervals, filename, trackname):
    """
    Exports a bed file track from dict with genomic positions.

    intervals=dict - chrom:(start, end) dict.
    filename=str.
    trackname=str.
    strand=str.
    """
    with open(filename, 'w') as handle:
        header = 'track name="{0}" description="{0}" visibility=pack autoScale=off colorByStrand="255,0,0 0,0,255"\n'.format(trackname)
        handle.write(header)
        for chrom, (start, end) in intervals.items():
            name = "DARN_{0}_{1}_{2}".format(chrom, start, end)
            entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                chrom, start, end,
                name,  # name
                1,  # score
                "."  # strand
            )
            handle.write(entry)


# Define variables
projectRoot = "/media/afr/cemm-backup/chipmentation"
data_dir = os.path.join(projectRoot, "periodicity")
results_dir = os.path.join(projectRoot, "periodicity")
plotsDir = os.path.join(results_dir, "/plots")

DNase = "/media/afr/cemm-backup/encode/wgEncodeUwDnaseK562Aln.merged.bam"
MNase = "/media/afr/cemm-backup/encode/wgEncodeSydhNsomeK562Aln.merged.bam"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "/chipmentation.replicates.annotation_sheet.csv"))

sampleSubset = samples[
    (samples["sampleName"].str.contains(
        "H3K4ME3_K562_500k_CM|" +
        "H3K4ME3_K562_500k_CHIP|" +
        "H3K27ME3_K562_500k_CM|" +
        "H3K27ME3_K562_500k_CHIP|" +
        "CTCF_K562_10M_CM|" +
        "PU1_K562_10M_CM|" +
        "CTCF_K562_10M_CHIP|" +
        "PU1_K562_10M_CHIP"
    ))
].reset_index(drop=True)

sampleSubset = sampleSubset.append(pd.Series(data=["DNase", DNase], index=["sampleName", "filePath"]), ignore_index=True)
sampleSubset = sampleSubset.append(pd.Series(data=["MNase", MNase], index=["sampleName", "filePath"]), ignore_index=True)

sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

# Temporary solution
bam_files = [
    "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam",
    "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam",
    "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam",
    "data/human/chipmentation/mapped/merged/H3K27me3_K562_10k500k_CM.bam",
    "data/human/chipmentation/mapped/merged/H3K27me3_K562_500k_ChIP.bam",
    "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam",
    "data/human/chipmentation/mapped/merged/CTCF_K562_10mio_CM.bam"
]

bedFilePath = "/home/arendeiro/reference/Homo_sapiens/hg19.cage_peak_coord_robust.TATA_Annotated.bed"
genome = "hg19"
windowRange = (-60, 60)
fragmentsize = 1
duplicates = True
n_clusters = 5

# plotly.sign_in("afrendeiro", "iixmygxac1")
windowWidth = abs(windowRange[0]) + abs(windowRange[1])

duplicates = False
window_width = 1000
window_step = 900
fragment_size = 1
genome = "hg19"


# Skip creating regions if exist (dedicated to Nathan)
regionsPickle = os.path.join(data_dir, "genomic_regions.no_repeats.pickle")
if os.path.isfile(regionsPickle):
    regions = pickle.load(open(regionsPickle, "r"))
else:
    # Get regions of interest in the genome
    gapsRepeats = BedTool(os.path.join("/home", "arendeiro", "reference/Homo_sapiens/hg19_gapsRepeats.bed"))

    # Whole genome in 1kb-windows
    whole_genome = BedTool.window_maker(BedTool(), genome='hg19', w=window_width, s=window_width)
    whole_genome = whole_genome.intersect(b=gapsRepeats, v=True, wa=True)

    # DHS and non-DHS in 1kb-windows
    dhs = BedTool(os.path.join("/home", "arendeiro", "wgEncodeOpenChromDnaseK562Pk.narrowPeak"))
    dhs = dhs.intersect(b=gapsRepeats, v=True, wa=True)
    non_dhs = whole_genome.intersect(b=dhs, v=True)

    # Bed files
    H3K27me3 = BedTool(os.path.join("/home", "arendeiro", "wgEncodeSydhHistoneK562H3k27me3bUcdPk.narrowPeak"))
    H3K27me3 = H3K27me3.intersect(b=gapsRepeats, v=True, wa=True)
    H3K4me3 = BedTool(os.path.join("/home", "arendeiro", "wgEncodeSydhHistoneK562H3k4me3bUcdPk.narrowPeak"))
    H3K4me3 = H3K4me3.intersect(b=gapsRepeats, v=True, wa=True)

    # peaks not overlapping the other mark in 1kb-windows
    H3K27me3_only = H3K27me3.intersect(b=H3K4me3, v=True)
    H3K4me3_only = H3K4me3.intersect(b=H3K27me3, v=True)

    # peaks overlapping each other in 1kb-windows
    H3K27me3_H3K4me3 = H3K4me3.intersect(b=H3K27me3)

    # Make 1kb windows and convert to HTSeq
    whole_genome = bedToolsInterval2GenomicInterval(whole_genome, strand=False, name=False)
    non_dhs = bedToolsInterval2GenomicInterval(BedTool.window_maker(non_dhs, b=non_dhs, w=window_width, s=window_width), strand=False, name=False)
    dhs = bedToolsInterval2GenomicInterval(BedTool.window_maker(dhs, b=dhs, w=window_width, s=window_width), strand=False, name=False)
    H3K27me3_only = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3_only, b=H3K27me3_only, w=window_width, s=window_width), strand=False, name=False)
    H3K4me3_only = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K4me3_only, b=H3K4me3_only, w=window_width, s=window_width), strand=False, name=False)
    H3K27me3_H3K4me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3_H3K4me3, b=H3K27me3_H3K4me3, w=window_width, s=window_width), strand=False, name=False)
    H3K27me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3, b=H3K27me3, w=window_width, s=window_width), strand=False, name=False)
    H3K4me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K4me3, b=H3K4me3, w=window_width, s=window_width), strand=False, name=False)

    regions = {"whole_genome": whole_genome,
               "dhs": dhs,
               "non_dhs": non_dhs,
               "H3K27me3": H3K27me3,
               "H3K4me3": H3K4me3,
               "H3K27me3_only": H3K27me3_only,
               "H3K4me3_only": H3K4me3_only,
               "H3K27me3_H3K4me3": H3K27me3_H3K4me3
               }
    pickle.dump(regions, open(os.path.join(data_dir, "genomic_regions.no_repeats.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

assert all([len(region) > 1 for region in regions.values()])

# with repeats, load:
# regions = pickle.load(open(os.path.join(data_dir, "genomic_regions.pickle"), "r"))

# Initialize Slurm object
slurm = DivideAndSlurm()
tasks = dict()

# Get sample names
samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in bam_files}

# Submit tasks for combinations of regions and bam files
for regionName, region in regions.items():
    for sampleName, sampleFile in samples.items():
        exportName = os.path.join(data_dir, sampleName + "_" + regionName)
        if os.path.isfile(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle")):
            continue
        # Add new task
        task = CountDistances(region, 4, os.path.abspath(sampleFile), permute=False, queue="shortq", time="4:00:00", permissive=True, cpusPerTask=4)
        slurm.add_task(task)
        slurm.submit(task)  # Submit new task
        tasks[task] = (sampleName, regionName, False)  # Keep track
        print(tasks[task])
        if os.path.isfile(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle")):
            continue
        # Add permuted
        task = CountDistances(region, 4, os.path.abspath(sampleFile), permute=True, queue="shortq", time="4:00:00", permissive=True, cpusPerTask=4)
        slurm.add_task(task)
        slurm.submit(task)  # Submit new task
        tasks[task] = (sampleName, regionName, True)  # Keep track
        print(tasks[task])

stored = list()
pickle.dump((slurm, tasks, stored), open("/home/arendeiro/slurm_20150210_1kb_norepeats.pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
slurm, tasks, stored = pickle.load(open("/home/arendeiro/slurm_20150210_1kb_norepeats.pickle", "r"))

# Collect processed data
for task, (sampleName, regionName, permuted) in tasks.items():          # loop through tasks, see if ready
    if task.is_ready():                                                 # if yes, collect output and save
        print(textwrap.dedent("""\
        Task {0} is now ready! {1}, {2}, {3}
        Time to completion was: {4} minutes.
        """.format(task, sampleName, regionName, permuted, int(time.time() - task.submissiontime) / 60.)))
        exportName = os.path.join(data_dir, sampleName + "_" + regionName)
        dists = task.collect()
        if not permuted:
            pickle.dump(dists, open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        else:
            pickle.dump(dists, open(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        stored.append(task)

# For each signal extract most abundant periodic signal through FFT, IFFT
for regionName, region in regions.items():
    for sampleName, sampleFile in samples.items():
        exportName = os.path.join(data_dir, sampleName + "_" + regionName)

        try:
            dists = pickle.load(open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "r"))
        except IOError("Can't open file."):
            continue

        distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
        distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}

        # Extract most abundant periodic pattern from signal
        # for DNase, extract from a different window (70-150bp)
        # patternPos = extractPattern(distsPos, range(60, 100), os.path.join(data_dir, exportName + "_posStrand-noRepeats"))
        # patternNeg = extractPattern(distsNeg, range(60, 100), os.path.join(data_dir, exportName + "_negStrand-noRepeats"))
        pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100), os.path.join(data_dir, exportName + "_bothStrands-noRepeats"))

        try:
            permutedDists = pickle.load(open(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle"), "r"))
        except IOError("Can't open file."):
            continue

        permutedDistsPos = {dist: count for dist, count in permutedDists.items() if dist >= 0}
        permutedDistsNeg = {abs(dist): count for dist, count in permutedDists.items() if dist <= 0}

        # Extract most abundant periodic pattern from signal
        # for DNase, extract from a different window (70-150bp)
        permutedPatternPos = extractPattern(permutedDistsPos, range(60, 100), os.path.join(data_dir, exportName + "_posStrand_permuted-noRepeats"))
        permutedPatternNeg = extractPattern(permutedDistsNeg, range(60, 100), os.path.join(data_dir, exportName + "_negStrand_permuted-noRepeats"))
        permutedPattern = extractPattern(Counter(permutedDistsPos) + Counter(permutedDistsNeg), range(60, 100), os.path.join(data_dir, exportName + "_bothStrands_permuted-noRepeats"))

# Focus on H3K4me3 data and nucleosome positioning
#
#
# calculate read coverage in H3K4me3 peaks
samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in bam_files}
sampleName = "H3K4me3_K562_500k_CM"
sampleFile = samples[sampleName]
regionName = "H3K4me3_only"
exportName = os.path.join(data_dir, sampleName + "_" + regionName)

# Correlate coverage and signal pattern
# get pattern
dists = pickle.load(open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "r"))
distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}
pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100),
                         os.path.join(data_dir, exportName + "_bothStrands-noRepeats"))

# Task get coverage and correlate with pattern, separate strands
slurm = DivideAndSlurm()
# make windows genome-wide with overlapping size of pattern
width = 1000
step = width - len(pattern)
genome_windows = BedTool.window_maker(BedTool(), g={"chr1": (1, 249250621)}, w=width, s=step)  # genome_windows = BedTool.window_maker(BedTool(), genome='hg19', w=width, s=step)
genome_windows = bedToolsInterval2GenomicInterval(genome_windows, strand=False, name=False)

task = CorrelatePatternBam(genome_windows.items(), 40, pattern, os.path.abspath(sampleFile), cpusPerTask=8)
slurm.add_task(task)
slurm.submit(task)
taskP = CorrelatePatternBam(genome_windows.items(), 40, pattern, os.path.abspath(sampleFile), permute=True, cpusPerTask=8)
slurm.add_task(taskP)
slurm.submit(taskP)

# collect and serialize
binary = task.collect()
assert len(genome_windows) == len(binary)  # compare sizes
assert all([len(window) == width - len(pattern) + 1 for window in binary.values()])  # likely to fail - patched downstream
pickle.dump(binary, open(os.path.join(data_dir, exportName + ".peakCorrelationBinary.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

binaryP = taskP.collect()
assert len(genome_windows) == len(binaryP)  # compare sizes
assert all([len(window) == width - len(pattern) + 1 for window in binaryP.values()])  # likely to fail - patched downstream
pickle.dump(binaryP, open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

binary = pickle.load(open(os.path.join(data_dir, exportName + ".peakCorrelationBinary.pickle"), "r"))
binaryP = pickle.load(open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryPermuted.pickle"), "r"))

# to subset, pass: OrderedDict(sorted((binary.items()[9000:10000]))) <- 1Mb of sequence
genome_binary = concatenateBinary(binary, len(pattern))  # 40
pickle.dump(genome_binary, open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryConcatenated.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

genome_binaryP = concatenateBinary(binaryP, len(pattern))
pickle.dump(genome_binaryP, open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryConcatenatedPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

genome_binary = pickle.load(open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryConcatenated.pickle"), "r"))
genome_binaryP = pickle.load(open(os.path.join(data_dir, exportName + ".peakCorrelationBinaryConcatenatedPermuted.pickle"), "r"))

# HMM
model = WinterHMM()

# Model training
#
# Train on subset of data, see probabilities
# i = len(genome_binary)
# model.train([genome_binary.values()[0][:10000]])  # subset data for training
# # see new probabilities
# [(s.name, s.distribution) for s in model.model.states]  # emission
# print(model.model.dense_transition_matrix())  # transition

# # save model with trained parameters
# pickle.dump(model, open(os.path.join(results_dir, sampleName + "_hmModel_trained_%i.pickle" % i), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
# model = pickle.load(open(os.path.join(results_dir, sampleName + "_hmModel_trained_%i.pickle" % i), "r"))

# Predict and get DARNS
hmmOutput = dict()
DARNS = dict()

# do it chromosome-wise
for chrom in genome_binary.keys():
    ite = range(0, len(genome_binary[chrom]), 5000000)
    prev = 0
    last = ite[-1]
    # do it for 5M chunks at a time
    for cur in ite[1:]:
        print(cur)
        if not cur == last:
            if chrom not in hmmOutput.keys():
                hmmOutput[chrom] = model.predict(genome_binary[chrom][prev:cur])
            else:
                hmmOutput[chrom] += model.predict(genome_binary[chrom][prev:cur])
        else:
            if chrom not in hmmOutput.keys():
                hmmOutput[chrom] = model.predict(genome_binary[chrom][prev:len(genome_binary[chrom])])
            else:
                hmmOutput[chrom] += model.predict(genome_binary[chrom][prev:len(genome_binary[chrom])])
        prev = cur

# add darns to dict
for chrom, sequence in hmmOutput.items():
    DARNS[chrom] = getDARNS(sequence)

pickle.dump((hmmOutput, DARNS), open(os.path.join(results_dir, sampleName + "_HMMResult.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
hmmOutput, DARNS = pickle.load(open(os.path.join(results_dir, sampleName + "_HMMResult.pickle"), "r"))

# Predict and get DARNS
hmmOutputP = dict()
DARNSP = dict()

# do it chromosome-wise
for chrom in genome_binaryP.keys():
    ite = range(0, len(genome_binaryP[chrom]), 5000000)
    prev = 0
    last = ite[-1]
    # do it for 5Mb chunks at a time
    for cur in ite[1:]:
        if not cur == last:
            if chrom not in hmmOutputP.keys():
                hmmOutputP[chrom] = model.predict(genome_binaryP[chrom][prev:cur])
            else:
                hmmOutputP[chrom] += model.predict(genome_binaryP[chrom][prev:cur])
        else:
            if chrom not in hmmOutputP.keys():
                hmmOutputP[chrom] = model.predict(genome_binaryP[chrom][prev:len(genome_binaryP[chrom])])
            else:
                hmmOutputP[chrom] += model.predict(genome_binaryP[chrom][prev:len(genome_binaryP[chrom])])
        prev = cur

# add darns to dict
for chrom, sequence in hmmOutputP.items():
    DARNSP[chrom] = getDARNS(sequence)

pickle.dump((hmmOutputP, DARNSP), open(os.path.join(results_dir, sampleName + "_HMMResultPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
hmmOutputP, DARNSP = pickle.load(open(os.path.join(results_dir, sampleName + "_HMMResultPermuted.pickle"), "r"))

# Filter DARNS with size = 1
DARNS = {chrom: [tup for tup in DARNS[chrom] if tup[1] - tup[0] > 1] for chrom in DARNS.keys()}
DARNSP = {chrom: [tup for tup in DARNSP[chrom] if tup[1] - tup[0] > 1] for chrom in DARNSP.keys()}

# Measure attributes:
# value of maximum positive strand correlation peak
# value of maximum negative strand correlation peak
# mean of positive strand correlation peak values
# mean of negative strand correlation peak values
# sum of positive strand correlation peak values
# sum of negative strand correlation peak values

# number of cycles through HMM “nucleosome” states (0.5 for each peak state passed)
# the HMM posterior probability of the path through the background state
# and the HMM posterior probability of the path through the nucleosome set of states

# length (bp)
# sum of reads within the DARNS
# read density
# spacing to nearest upstream DARNS (bp)
# spacing to nearest downstream DARNS (bp)

# maybe also read count and density from MNase

# Measure read count and density in DARNS

# in serial
cov = dict()
for regions in ["DARNS", "DARNSP"]:
    intervals = DARNS2GenomicInterval(eval(regions))
    bam = HTSeq.BAM_Reader("/media/afr/cemm-backup/chipmentation/data/mapped/K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam")
    cov[regions] = coverage(
        bam,
        intervals,
        1
    )

# in parallel
slurm = DivideAndSlurm()

for regions in ["DARNS", "DARNSP"]:
    # make windows genome-wide with overlapping size of pattern
    intervals = DARNS2GenomicInterval(eval(regions))

    task = Coverage(intervals, 40, os.path.abspath(samples["H3K4me3_K562_500k_CM"]), cpusPerTask=8)
    task.id = regions
    slurm.add_task(task)
    slurm.submit(task)

while not all([t.is_ready() for t in slurm.tasks]):
    time.sleep(10)

# collect and serialize
cov = dict()
for task in slurm.tasks:
    cov[task.id] = task.collect()

# serialize coverage
pickle.dump(
    cov,
    open(os.path.join(results_dir, sampleName + "_DARNS.coverage.pickle"), "wb"),
    protocol=pickle.HIGHEST_PROTOCOL
)

cov = pickle.load(open(os.path.join(results_dir, sampleName + "_DARNS.coverage.pickle"), 'r'))

# Measure DARNS features in parallel
# using tasks across nodes
features = ["length", "space_upstream", "space_downstream",
            "read_count", "read_density", "n_pospeaks", "n_negpeaks"]  # anyway already hard-coded in task

mode = "parallel"

if mode == "serial":
    # in serial
    import functools
    DARNS_features = pd.DataFrame()

    for region in ['DARNS', "DARNSP"]:
        for chrom in eval(region).keys():
            df = pd.DataFrame(  # this is the reduce step
                map(
                    functools.partial(  # use wrapper to pass more arguments
                        measureDarnsFeatures,
                        eval(region),
                        region,
                        chrom,
                        region,
                        features,
                        cov
                    ),
                    range(len(eval(region)[chrom]))
                )
            )
            DARNS_features = DARNS_features.append(df, ignore_index=True)
elif mode == "parallel":
    # in parallel
    import multiprocessing
    import parmap
    DARNS_features = pd.DataFrame()

    for region in ['DARNS', "DARNSP"]:
        for chrom in eval(region).keys():
            df = pd.DataFrame(  # this is the reduce step
                parmap.map(
                    measureDarnsFeatures,
                    range(len(eval(region)[chrom])),  # iterable
                    regions=eval(region),
                    hmm=region,
                    chrom=chrom,
                    regionType=region,
                    features=features,
                    coverage=cov
                )
            )
            DARNS_features = DARNS_features.append(df, ignore_index=True)
elif mode == "slurm":
    # in parallel across nodes
    slurm = DivideAndSlurm()

    for region in ['DARNS', "DARNSP"]:
        for chrom in eval(region).keys():
            task = MeasureDarnsFeatures(
                range(len(eval(region))),  # indexes of darns
                100,  # n of jobs to split across
                chrom,
                region
            )
            slurm.add_task(task)
            slurm.submit(task)

    DARNS_features = pd.DataFrame()
    for task in slurm.tasks:
        df = task.collect()
        # append
        DARNS_features = DARNS_features.append(df, ignore_index=True)

# Remove NAs (only first and last darn)
DARNS_features = DARNS_features.dropna().reset_index(drop=True)

# Fill in data by imputation
# from sklearn.preprocessing import Imputer
# imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
# imp.fit(DARNS_features.dropna())  # need to subset and transform to np.array
# DARNS_features = imp.transform(DARNS_features)

# serialize
pickle.dump(
    DARNS_features,
    open(os.path.join(results_dir, sampleName + "_DARNS.features.pickle"), "wb"),
    protocol=pickle.HIGHEST_PROTOCOL
)
# DARNS_features = pickle.load(open(os.path.join(results_dir, sampleName + "_DARNS.features.pickle"), "r"))
DARNS_features = pd.read_pickle(os.path.join(results_dir, sampleName + "_DARNS.features.pickle"))

# Test significant differences between Real and Permuted DARNS
from scipy.stats import ks_2samp
from scipy.stats import chisquare
for feature in features:
    real = DARNS_features.loc[DARNS_features["type"] == "DARNS", feature].tolist()
    permuted = DARNS_features.loc[DARNS_features["type"] == "DARNSP", feature].tolist()

    # test difference
    D, p = ks_2samp(real, permuted)
    print(feature, D, p)
    # see normality
    # qqplot
    import scipy.stats as stats

    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    f.suptitle('%s - p:%f' % (feature, p))

    stats.probplot(real, dist="norm", plot=ax1)
    stats.probplot(permuted, dist="norm", plot=ax2)

    ax1.set_title('%s - Real data' % feature)
    ax2.set_title('%s - Permuted data' % feature)

    # f.tight_layout()
    plt.savefig(os.path.join(results_dir, sampleName + "_DARNS.qqplot-%s.pdf" % feature))
    plt.close()

    # plot variable distribution
    sns.set(style="white", palette="muted")
    b, g, r, p = sns.color_palette("muted", 4)
    f, axes = plt.subplots(4, 2, figsize=(7, 7), sharex=True)

    sns.distplot(real, kde=False, color=b, ax=axes[0, 0])
    sns.distplot(real, hist=False, rug=True, color=r, ax=axes[0, 1])
    sns.distplot(real, hist=False, color=g, kde_kws={"shade": True}, ax=axes[1, 0])
    sns.distplot(real, color=p, ax=axes[1, 1])

    sns.distplot(permuted, kde=False, color=b, ax=axes[2, 0])
    sns.distplot(permuted, hist=False, rug=True, color=r, ax=axes[2, 1])
    sns.distplot(permuted, hist=False, color=g, kde_kws={"shade": True}, ax=axes[3, 0])
    sns.distplot(permuted, color=p, ax=axes[3, 1])

    plt.savefig(os.path.join(results_dir, sampleName + "_DARNS.distplot-%s.pdf" % feature))
    plt.close()

# plt.show()

# Plot features for Real and Permuted DARNS
# corrplot between all variables
sns.set(style="darkgrid")

f, ax = plt.subplots(figsize=(9, 9))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.corrplot(DARNS_features[features], annot=False, sig_stars=False,
             diag_names=False, cmap=cmap, ax=ax)
f.tight_layout()
plt.savefig(os.path.join(results_dir, sampleName + "_DARNS.corrplot.pdf"))
plt.close()

# pair plot
f = sns.PairGrid(DARNS_features[features], diag_sharey=False)
f.map_lower(sns.kdeplot, cmap="Blues_d")
f.map_upper(plt.scatter)
f.map_diag(sns.kdeplot, lw=3)
plt.savefig(os.path.join(results_dir, sampleName + "_DARNS.pairplot.pdf"))
plt.close()

# pairwise distributions
size = sum(1 for _ in itertools.combinations(features, 2))
f, axs = plt.subplots(size, sharex=True, sharey=True)
i = 0
for f1, f2 in itertools.combinations(features, 2):
    axs[i] = sns.jointplot(DARNS_features[f1], DARNS_features[f2], kind="kde", size=7, space=0)
    i += 1
plt.savefig(os.path.join(results_dir, sampleName + "_DARNS.pairplot.pdf"))
plt.close()

# Train classifier on features
# You can start from here.
# Previous variables:
samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in bam_files}
sampleName = "H3K4me3_K562_500k_CM"
sampleFile = samples[sampleName]
regionName = "H3K4me3_only"
exportName = os.path.join(data_dir, sampleName + "_" + regionName)
features = ["length", "space_upstream", "space_downstream",
            "read_count", "read_density", "n_pospeaks", "n_negpeaks"]
features = ["length", "space_upstream", "space_downstream",
            "read_count", "read_density", "n_pospeaks", "n_negpeaks"]

DARNS_features = pd.read_pickle(os.path.join(results_dir, sampleName + "_DARNS.features.pickle"))

# Standardize
# from sklearn import preprocessing
# for feature in features:
#     DARNS_features.loc[:, feature] = preprocessing.scale(DARNS_features[feature])

# get random tenth of data to train on
randomRows = [random.randrange(0, len(DARNS_features)) for _ in range(len(DARNS_features) / 10)]
train_features = DARNS_features.loc[randomRows, features]

# get class labels for data
train_labels = DARNS_features.loc[randomRows, "type"].tolist()

for classifier in ["linear", "tree", "forest", "neighbours"]:
    print(classifier)
    if classifier == "linear":
        # Linear model classifier - logistic regressor
        from sklearn.linear_model import LogisticRegression
        clf = LogisticRegression()
    elif classifier == "svm":
        from sklearn import svm
        clf = svm.SVC()
    elif classifier == "neighbours":
        # Use nearest neighbours classifier
        from sklearn import neighbors
        n_neighbors = 4
        clf = neighbors.KNeighborsClassifier(n_neighbors, weights="uniform")  # try with "distance"
    elif classifier == "tree":
        # Use decision tree
        from sklearn.tree import DecisionTreeClassifier
        clf = DecisionTreeClassifier()
    elif classifier == "forest":
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier()

    # Train
    print("training...")
    clf.fit(train_features, train_labels)

    # Predict for all data
    # put in data frame form
    # pred = pd.DataFrame()
    # pred["name"] = DARNS_features["name"].tolist()
    # pred["type"] = DARNS_features["type"].tolist()
    # split = list(DARNS_features["name"].apply(lambda x: x.split("_")))
    # pred["chr"] = [i[0] for i in split]
    # pred["start"] = [i[1] for i in split]
    # pred["end"] = [i[2] for i in split]
    # pred["center"] = [i[3] for i in split]

    # Predict N times for all data
    N = 10

    pred = list()
    for i in range(N):  # try increasing N
        print("predicting ", i)
        pred += list(clf.predict_proba(DARNS_features.loc[:, features])[:, 0])

    # pickle.dump(
    #     pred,
    #     open(os.path.join(results_dir, sampleName + "_DARNS." + classifier + ".predictions.pickle"), "wb"),
    #     protocol=pickle.HIGHEST_PROTOCOL
    # )

    # Get metrics
    print("getting metrics")
    y_true = [0 if i == "DARNSP" else 1 for i in DARNS_features["type"].tolist()] * N

    # report
    # classification_report(y_true, pred)

    # ROC
    fpr, tpr, thresholds = roc_curve(y_true, pred)
    AUC = roc_auc_score(y_true, pred)

    # PRC
    precision, recall, thresholds = precision_recall_curve(y_true, pred)
    APS = average_precision_score(y_true, pred)

    print("plotting")
    # plot ROC
    plt.plot(fpr, tpr)
    plt.text(0.05, 0.95, "AUC = %s" % AUC, fontsize=12)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.savefig(os.path.join(results_dir, sampleName + "_DARNS." + classifier + ".roc.pdf"))
    plt.close()

    # plot PRC
    plt.plot(recall, precision)
    plt.text(0.05, 0.95, "AUC = %s" % APS, fontsize=12)
    plt.plot([0, 1], [1, 0], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.savefig(os.path.join(results_dir, sampleName + "_DARNS." + classifier + ".prc.pdf"))
    plt.close()


# Calculate FDR for each DARN

# Export DARN tracks
allDARNS = pred[pred["type"] == "DARNS"]

# invert FDR and scale 0-1000 for bed track
allDARNS.replace({"FDR": {np.inf: 1}}, inplace=True)
allDARNS.loc[:, "FDR"] = (1 - allDARNS["FDR"]) * 1000

allDARNS.loc[:, ["chr", "start", "end", "name", "FDR"]].to_csv(
    os.path.join(results_dir, sampleName + "_DARNS.all.bed"),
    sep="\t", index=False, header=False
)

# Keep DARNS with FDR < 0.5
len(pred[pred["FDR"] < 0.5])
realOnes = pred[(pred["FDR"] < 0.5) & (pred["type"] == "DARNS")]

# invert FDR and scale 0-1000 for bed track
realOnes.replace({"FDR": {np.inf: 1}}, inplace=True)
realOnes.loc[:, "FDR"] = (1 - realOnes["FDR"]) * 1000

realOnes.loc[:, ["chr", "start", "end", "name", "FDR"]].to_csv(
    os.path.join(results_dir, sampleName + "_DARNS.FDR_filtered.bed"),
    sep="\t", index=False, header=False
)

fakeOnes = pred[pred["type"] == "DARNSP"]
# invert FDR and scale 0-1000 for bed track
fakeOnes.replace({"FDR": {np.inf: 1}}, inplace=True)
fakeOnes.loc[:, "FDR"] = (1 - fakeOnes["FDR"]) * 1000

fakeOnes.loc[:, ["chr", "start", "end", "name", "FDR"]].to_csv(
    os.path.join(results_dir, sampleName + "_DARNSP.all.bed"),
    sep="\t", index=False, header=False
)

# Plot again variables for DARNS with FDR < 0.5

# Alternatively, train only with real DARNS overlapping H3K4me3 peaks

# Plot:
#    - From the middle of the DARN (mid-peak)
#    - From the 5' and 3' end of nucleosome Dyads
#    - DARNS frequency around TSSs (models and CAGE) and TTSs
#    - DARNS frequency around CpGs islands





##### OLD #####

# Plot attributes
assert len(DARNS) == len(DARNSP)
for data in (DARNS, DARNSP):
    widths = Counter([darn[1] - darn[0] for darn in data.values()[0]])
    distances, midDistances = measureDARNS(data)

    plt.plot(widths.keys(), widths.values(), '-', color='orange')
    plt.plot(distances.keys(), distances.values(), '-', color='purple')
    plt.plot(midDistances.keys(), midDistances.values(), '-', color='green')
    # sns.violinplot(widths.values())
    # decide how to save

# Get scores from DARNS
probs = {peak: model.retrieveProbabilities(sequence) for peak, sequence in genome_binary.items()}

# Plot predicted darns and post probs in an example region
#
chrom = "chr1"
start = 5990
end = 6000
width = (end - start) * 1000  # length in bp

# retrieve post. probs for selected region
probs = model.retrieveProbabilities(OrderedDict(sorted((genome_binary[chrom][start:end]))))  # random 10kb region
probsP = model.retrieveProbabilities(OrderedDict(sorted((genome_binaryP[chrom][start:end]))))  # random 10kb region

# start plotting
from matplotlib.patches import Rectangle
colors = sns.color_palette('deep', n_colors=6, desat=0.5)
sns.set_context(rc={"figure.figsize": (14, 6)})
sns.plt.axhline(y=1.1, c=colors[0], alpha=0.7)
sns.plt.xlim([1, width + 1])
sns.plt.ylim([0, 1.2])
sns.plt.ylabel(r'posterior probs, $\gamma_k$')
sns.plt.xlabel(r'$k$')
axis = sns.plt.gca()

# viterbi predicted DARNS
# get darns in window
DARNSinWindow = [darn for darn in DARNS[chrom] if darn[0] >= start and darn[1] <= end]
# for each darn in window draw a box
for start, end in DARNSinWindow:
    axis.add_patch(Rectangle((start + 1, 1.075), end - start + 1, 0.05,
                             facecolor=colors[0], alpha=0.7))

# line plot of post. probs
sns.plt.plot(range(1, width + 1), probs,  # post. probs of real data
             c=colors[2], alpha=0.7)
sns.plt.plot(range(1, width + 1), probsP,  # post. probs of permuted data
             c=colors[3], alpha=0.5)
plt.show()  # decide how to save

# TO IMPLEMENT:
# Export bed files of predicted DARNS for both real data and permuted
exportBedFile(
    DARNS,
    os.path.join(results_dir, sampleName + ".DARNS.bed"),
    "DARNS predicted from %s" % sampleName
)
exportBedFile(
    DARNS,
    os.path.join(results_dir, sampleName + ".DARNS.permuted.bed"),
    "DARNS predicted from %s permuted reads" % sampleName
)
# Get overal score for darn from post. prob (or something else)

# Get overal score for darn over permuted

# Output bed/wig with scores!
# Export wig files with raw correlations
# exportWigFile(
#     [peaks[i] for i in correlationsPos.keys()],
#     correlationsPos.values(),
#     len(pattern) / 2,
#     os.path.join(results_dir, sampleName + ".peakCorrelationPos.wig"),
#     sampleName + " raw absolute correlation - positive strand"
# )
# exportWigFile(
#     [peaks[i] for i in correlationsNeg.keys()],
#     correlationsNeg.values(),
#     len(pattern) / 2,
#     os.path.join(results_dir, sampleName + ".peakCorrelationNeg.wig"),
#     sampleName + " raw absolute correlation - negative strand"
# )

# TODO:
# Get regions in extreme quantiles of correlation
# Check for enrichment in ...
#     mnase signal
#     nucleosomes
#     clusters of CM signal around TSSs
