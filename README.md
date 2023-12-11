# TX2P
TX2P - Transcript to protein


## Getting Started:
Download the docker images.

>docker pull murphydaviducl/getorf:latest
>
>docker pull murphydaviducl/metamorpheusdocker:latest
>
>git clone https://github.com/MurphyDavid/TX2P
>
>cd TX2P

"metamorpheusdocker" is a copy of the one availible from https://smith-chem-wisc.github.io/MetaMorpheus/ or  https://github.com/smith-chem-wisc/MetaMorpheus
If desired you should be able to use an updated version from those sources. Please cite https://doi.org/10.1021/acs.jproteome.7b00873


## Testing:

>current_folder=$(pwd)

>docker run -v $current_folder:/current_folder  --entrypoint Rscript murphydaviducl/getorf /getcds/getCDSandAAseq.R --gtf /current_folder/test.gtf --transcript /current_folder/test_transcripts.tsv --output /current_folder/output_ORFs.csv


