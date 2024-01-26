Download and install docker

 
![image](https://github.com/MurphyDavid/TX2P/assets/11276387/b33b53fd-f93f-445b-bc3c-33690cbfebdb)



Once it’s installed, create a docker account and then search for images to run

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/ac7cd367-786d-46d7-96c0-e9a3d97ac40a)

 
Click pull for murphydaviducl/getorf:latest

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/9a1b94d3-a1e2-4bc4-bd91-570decca18fb)

Click pull for murphydaviducl/metamorpheusdocker:latest
 
![image](https://github.com/MurphyDavid/TX2P/assets/11276387/969c5fb8-34ee-43d3-b58d-7277e9409fe0)

Download the TX2P code from github : https://github.com/MurphyDavid/TX2P/

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/4a6ddb75-4cd6-4230-ab71-2556877fe634)

Extract the zip file to a folder.

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/28b9ae07-6a24-40d1-98a2-c3d95990f08c)

 
Run CMD

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/3a2a975c-3dd3-440c-aba7-b27ddae31380)
![image](https://github.com/MurphyDavid/TX2P/assets/11276387/a483f460-9ffb-476a-846d-2a27f50281df)


Go to the folder in CMD

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/db46307d-5791-4544-a876-3f40976af22f)


Download the example mass spec raw file with this curl command:

```
curl ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020044/H_Luh_ND_3.raw
 --output H_Luh_ND_3.raw
```

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/e14f98c7-4536-48e4-9ae0-346513582db3)


```
set current_folder="%cd%"
echo %current_folder%
```

Now run the get-orf tool

```
docker run -v %current_folder%:/current_folder  --entrypoint Rscript murphydaviducl/getorf /getcds/getCDSandAAseq.R --gtf /current_folder/test.gtf --transcript /current_folder/test_transcripts.tsv --output /current_folder/output_ORFs.csv
```

 ![image](https://github.com/MurphyDavid/TX2P/assets/11276387/154392a1-ec82-42d2-a840-ff93d4203467)

```
mkdir testoutput
cd ./testoutput

docker run -v %current_folder%:/current_folder -it --entrypoint dotnet murphydaviducl/metamorpheusdocker /metamorpheus/CMD.dll -o /current_folder/testoutput -t /current_folder/mmconfig/Task2-CalibrateTaskconfig.toml /current_folder/mmconfig/Task4-GPTMDTaskconfig.toml current_folder/mmconfig/Task5-SearchTaskconfig.toml -d /current_folder/test.fasta -s /current_folder/H_Luh_ND_3.raw
```

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/169388c1-2df5-4512-863e-43b0ccbeebc3)

 
This example will take a little while to run, running with more seuquences and more mass spec data may take much longer.

Using the example files provided, once the tool has run you should now have a folder that looks like this.


![image](https://github.com/MurphyDavid/TX2P/assets/11276387/4b54524f-4654-4cbc-8124-6b04e7e956e0)

The Task3 folder should look like this

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/79156715-18f9-43cd-a1fb-e5e93057300c)

The AllQuantifiedProteinGroups should have one entry like this

![image](https://github.com/MurphyDavid/TX2P/assets/11276387/cb8e5840-9ca7-4aa0-aeeb-73bace0465fd)


