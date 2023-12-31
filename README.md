<p align="center">
  <img src="img/RNA_protein.png" width="500" height="400"/>  
</p>

# TX2P - Transcript to protein

<!-- badges: start -->
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`TX2P` allows for automated integration of mass spectometry data in long-read RNA-sequencing workflow. It does so by predicting open reading frames for each transcripts within an input GTF/GFF file, translates it into peptide sequences and then searches for those predicted proteins in mass spec datasets using `MetaMorpheus`.

*NOTE!* Currently `TX2P` is only set up to work with GRCh38.

## Prerequisites

Before installing `TX2P`, you need to have `Docker` installed on your system. `Docker` is a platform that allows you to create, deploy, and run applications in containers. Follow the instructions below to install `Docker` if you haven't already:

### Installing Docker

#### For Linux:

1. **Open a Terminal:** Use your Linux distribution's package manager to install `Docker`. The command varies depending on the distribution:
   **Ubuntu/Debian**: 
   ```
   sudo apt-get install docker-ce docker-ce-cli containerd.io
   ```
   **Fedora**: 
   ```
   sudo dnf -y install docker-ce
   ```
   **CentOS**: 
   ```
   sudo yum install docker-ce docker-ce-cli containerd.io
   ```
2. **Start the Docker Service:** Use `sudo systemctl start docker` to start the Docker service.
3. **Verify Installation:** Check if Docker is installed correctly by typing `docker --version` in the terminal.

#### For Windows

The docker's will work on windows but TODO: exact shell commands

## Getting Started:

Download the docker images.

```
docker pull murphydaviducl/getorf:latest

docker pull murphydaviducl/metamorpheusdocker:latest

git clone https://github.com/MurphyDavid/TX2P

cd TX2P
```

`metamorpheusdocker` is a copy of `MetaMorpheus` as availible from https://smith-chem-wisc.github.io/MetaMorpheus/ or  https://github.com/smith-chem-wisc/MetaMorpheus
If desired you should be able to use an updated version from those sources. 
Please cite https://doi.org/10.1021/acs.jproteome.7b00873


## Testing:

### Getting ORFs

```
current_folder=$(pwd)

docker run -v $current_folder:/current_folder  --entrypoint Rscript murphydaviducl/getorf /getcds/getCDSandAAseq.R --gtf /current_folder/test.gtf --transcript /current_folder/test_transcripts.tsv --output /current_folder/output_ORFs.csv
```

### Searching mass spec data.

#### Download a mass spec file for testing.

```
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020044/H_Luh_ND_3.raw
```
#### Create an output folder
```
current_folder=$(pwd)

mkdir testoutput

cd ./testoutput
```
#### Run on the example mass spec data. 

```
docker run -v $current_folder:/current_folder -it --entrypoint dotnet mymodifiedmetamorpheus:latest /metamorpheus/CMD.dll -o /current_folder/testoutput -t /current_folder/mmconfig/Task2-CalibrateTaskconfig.toml /current_folder/mmconfig/Task4-GPTMDTaskconfig.toml current_folder/mmconfig/Task5-SearchTaskconfig.toml -d /current_folder/test.fasta -s /current_folder/H_Luh_ND_3.raw
```
You can add almost any number of additional mass spec files, it's recommended to include a full dataset to allow good calibration though we suggest keeping it to less than 200GB or so. 
You can give the docker access to additional folders by adding  more in the format "-v folder:/mountname" 
The paths to additional raw files needs to be give in terms of the mounted path inside the docker container. 

Supported formats are 

>Thermo .raw 
>
>.mzML file in centroid mode.
>
>.mgf

Using the example files provided, once the tool has run you should now have a folder that looks like this.

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/c899d7b9-80d2-472a-b4bf-a6aae6b6b1ad)

The Task3 folder should look like this

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/c6a4eeab-2faf-42a2-b305-ca64983e2082)

The AllQuantifiedProteinGroups should have one entry like this

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/0717f0ad-189b-4911-bd6d-dc43a82cd755)

### Usage

If you are trying to confirm that a transcript in an organism is producing a protein, it's suggested to include the full known proteome for that organism just as you would a contaminants file.
For humans the relevant fasta file can be downloaded from uniprot: https://www.uniprot.org/help/downloads 

>docker run -v $current_folder:/current_folder -it --entrypoint dotnet mymodifiedmetamorpheus:latest /metamorpheus/CMD.dll -o /current_folder/testoutput -t /current_folder/mmconfig/Task2-CalibrateTaskconfig.toml /current_folder/mmconfig/Task4-GPTMDTaskconfig.toml current_folder/mmconfig/Task5-SearchTaskconfig.toml -d [Fasta file for transcripts of interest] /current_folder/contaminants.fasta [Fasta file of known proteins for organism of interest] -s [complete list of mass spec files]
