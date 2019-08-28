import sys
import os
import subprocess

def execute(command):

    print('Executing {}'.format(command))
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    return stdout, stderr

output_dir = snakemake.output

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

identity = snakemake.params.identity
coverage = snakemakem.params.coverage
cluster_cmd = './scripts/Cluster_genomes.pl -f {0} -c {1} -i {2}'.format(os.path.basename(snakemake.input), coverage, identity)

cluster_stdout, cluster_stderr = execute(cluster_cmd)

with open('ClusterGenomes.err', 'w') as stderr_fh:
    stderr_fh.write(cluster_stderr + os.sep)

with open('ClusterGenomes.out', 'w') as stdout_fh:
    stdout_fh.write(cluster_stdout + os.sep)

