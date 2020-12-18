# Download docker and pull the docker image of pasa from:
# https://github.com/PASApipeline/PASApipeline/wiki/PASA_Docker
# Make a new container and access through the command line (CLI button in docker)

#---------------------------------------------
# Files needed
# trinity de novo: trinity_annotations_sex_specific.fasta
# trinity genome guided: Trinity-GG.fasta
# annotation: PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3
# genome: PCITRI.assembly.v0.fa 

#---------------------------------------------
# Following this pipeline:
# https://github.com/PASApipeline/PASApipeline/wiki/PASA_comprehensive_db

#https://stackoverflow.com/questions/22907231/how-to-copy-files-from-host-to-docker-container
#https://docs.windriver.com/bundle/Wind_River_Linux_Tutorial_Using_Docker_Containers_CD/page/kri1507808414748.html

echo "copy in data to dockercontainer using normal terminal"
docker ps # to get container ID (change ID below to correct one)
docker cp trinity_annotations_sex_specific.fasta 416e49a0c862:/trinity_annotations_sex_specific.fasta
docker cp Trinity-GG.fasta 416e49a0c862:/Trinity-GG.fasta
docker cp PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3 416e49a0c862:/PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3
docker cp PCITRI.assembly.v0.fa 416e49a0c862:/PCITRI.assembly.v0.fa
# in docker CLI mv all files to /usr/local/src (directory where pasa is so everything will run)

echo "in docker CLI: make config files"
# inside the pasa_conf directory which comes with PASA
cp pasa.alignAssembly.Template.txt alignAssembly.config # change line to read: DATABASE=/tmp/my_database
cp pasa.annotationCompare.Template.txt annotCompare.config # change line to read: DATABASE=/tmp/my_database
# mv both files to the main container: /usr/local/src

echo "in docker CLI: format transcriptome annotations"
cat trinity_annotations_sex_specific.fasta Trinity-GG.fasta > transcripts.fasta
./PASApipeline/misc_utilities/accession_extractor.pl < trinity_annotations_sex_specific.fasta > tdn.accs

echo "in docker CLI: run pasa"
./PASApipeline/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R --TDN tdn.accs -g PCITRI.assembly.v0.fa --ALIGNERS blat -t transcripts.fasta --transcribed_is_aligned_orient -L --annots PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3 --gene_overlap 50.0

echo "in docker CLI: make transcriptome database"
./PASApipeline/scripts/build_comprehensive_transcriptome.dbi -c alignAssembly.config -t transcripts.fasta --min_per_ID 95 --min_per_aligned 30

echo "in docker CLI: make final file - round 1"
./PASApipeline/misc_utilities/pasa_gff3_validator.pl PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3
./PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g PCITRI.assembly.v0.fa -P PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3
./PASApipeline/Launch_PASA_pipeline.pl -c annotCompare.config -A -g PCITRI.assembly.v0.fa -t transcripts.fasta

echo "in docker CLI: make final file - round 2 using the previously made new .gff3"
./PASApipeline/misc_utilities/pasa_gff3_validator.pl my_database.gene_structures_post_PASA_updates.132.gff3
./PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g PCITRI.assembly.v0.fa -P my_database.gene_structures_post_PASA_updates.132.gff3
./PASApipeline/Launch_PASA_pipeline.pl -c annotCompare.config -A -g PCITRI.assembly.v0.fa -t transcripts.fasta

echo "in terminal: mv files to my computer"
docker cp 416e49a0c862:/usr/local/src/my_database.gene_structures_post_PASA_updates.15526.gff3 ./my_database.gene_structures_post_PASA_updates.15526.gff3
