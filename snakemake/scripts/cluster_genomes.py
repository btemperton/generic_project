#!/usr/bin/env python

###############################################################################
#                                                                             #
#    Cluster_genomes                                                          #
#                                                                             #
#    A wrapper script for Cluster_genomes.pl and written for Docker for use   #
#    with Cyverse's Docker platform.                                          #
#                                                                             #
#    Copyright (C) Benjamin Bolduc                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

__author__ = "Ben Bolduc"
__copyright__ = "Copyright 2016"
__credits__ = ["Ben Bolduc"]
__license__ = "LGPLv3"
__maintainer__ = "Ben Bolduc"
__email__ = "bolduc.10@osu.edu"
__status__ = "Development"

import sys
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(
    description='Clusters genomes',
    formatter_class=argparse.RawTextHelpFormatter)

options = parser.add_argument_group('Options')
options.add_argument('-f', '--fasta', dest='fasta_fn', metavar='DIRECTORY',
                     help='Filename of fasta file to cluster.')
options.add_argument('-c', '--coverage', dest='coverage', default=80,
                     help='Percent of the sequence covered.')
options.add_argument('-i', '--identity', dest='identity', default=95,
                     help='Percent identity')
options.add_argument('-o', '--output-dir', dest='output_dir', default='Output',
                     help='Directory to place output files')

results = parser.parse_args()


def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)


def execute(command):

    print('Executing {}'.format(command))
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    return stdout, stderr


if __name__ == '__main__':

    cwd = os.getcwd()
    output_dir = os.path.join(cwd, results.output_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cp_inputs = 'cp {0} {1}'.format(results.fasta_fn, output_dir)
    execute(cp_inputs)

    os.chdir(output_dir)

    identity = results.identity
    coverage = results.coverage
    cluster_cmd = 'Cluster_genomes.pl -f {0} -c {1} -i {2}'.format(os.path.basename(results.fasta_fn), coverage, identity)

    cluster_stdout, cluster_stderr = execute(cluster_cmd)

    with open('ClusterGenomes.err', 'w') as stderr_fh:
        stderr_fh.write(cluster_stderr + os.sep)

    with open('ClusterGenomes.out', 'w') as stdout_fh:
        stdout_fh.write(cluster_stdout + os.sep)
        
    try: 
        os.remove(os.path.join(cwd, os.basename(results.fasta_fn)))  # chdir earlier, so it's only the duplicate, not the original
    except Exception as e:
        print('Error with removing file {}. This isnt a problem if the output directory contains *.clstr file. If it doesnt, check the *.err log file for any errors.'.format(os.basename(results.fasta_fn)))
    print('Program Complete')

