rule initialise_virsorter:
  output:
    home=dir("tools/VirSorter"),
    data=dir("tools/VirSorter/virsorter-data")
  conda: "../envs/virsorter.yml"
  shell:
    """
      cd tools;
      git clone https://github.com/simroux/VirSorter.git;
      cd VirSorter/Scripts;
      make clean;
      make;
      cd ..;
      wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz;
      tar -xvzf virsorter-data-v2.tar.gz;
      rm virsorter-data-v2.tar.gz;
      cd ../../;
    """

rule initialise_DeepVirFinder:
  output: 
    cmd="tools/DeepVirFinder/dvf.py",
    directory="tools/DeepVirFinder"
  conda: "../envs/DeepVirFinder.yml"
  shell:
    """
    cd tools;
    git clone https://github.com/jessieren/DeepVirFinder;
    cd DeepVirFinder;
    chmod +x dvf.py;
    cd ../../;
    """
