# De novo transcriptome assembly

Here I make a de novo transcriptome assembly from the male and female RNA-Seq data in order to better annotate alternative splicing. (This was suggested by reviewers).

Following Andrés's pipeline found [here](https://github.com/RossLab/B_viburni/blob/master/2_Transcriptome.md) to make the trinity assmebly. 

---

Now I have an assembly I will try blasting this against the genome to get a the corresponding gene_ids so I can rename the trinity assembly transcripts with P.citri_v0 gene_ids, this will allow me to do all downstream analysis with the methylation data.

After a quick blast the intitial result gives the Pcitri genes in terms of the location so need to match up the location to the gene name (have this in a previous script though somewhere).

---

This doesn't work, found [this](https://github.com/PASApipeline/PASApipeline/wiki) though which 'should' do what I need, annotating alt splicing but keep in the appropriate gene ids. It's recommended by Trinity and needs two Trinity inputs, a de novo transcriptome and a genome guided de novo transcriptome.