# use these 2 commands to run the pipeline and post-process the data

# First:
# configure virtualBox to mount folders
# ./VBoxManage.exe sharedfolder add boot2docker-vm --name "d/" --hostpath "d:/" --automount
# ./VBoxManage.exe sharedfolder add boot2docker-vm --name "UDG/" --hostpath "M:\DB_Unite de Genomique\Poletti_NGS data" --automount
# ./VBoxManage.exe sharedfolder add boot2docker-vm --name "DBIMBI/" --hostpath "M:\DB_IMBI\Bioinformatique" --automount


# then map local folder (here it is d:/) to boot2docker:
# sudo mount -t vboxsf d/ /home/
# or 
# sudo mount -t vboxsf DBIMBI/ /home/


# then :
# run these commands to run the pipeline and any post-processing scripts like annotation or reporting.

docker run --rm -v /home/:/home/ -w $(pwd)  genethon/fragvisa sh /home/tempGC/Scripts/start_snake.sh $PWD;
docker run --rm -v /home/:/home/ -w $(pwd)  gc/stat_is Rscript /home/tempGC/Scripts/Master_r.r $(pwd);

# Explanations:
# run the docker image "fragvisa", with folder mapping home(boot2docker which was mapped to d:) to the home folder of the container. Every path is now transparent.
# Set the working directory in the container to the current directory in boot2docker.
# Run the shell script that contains snakemake commands.
#
# Second command is the same but for Rscript instead of snakemake
#