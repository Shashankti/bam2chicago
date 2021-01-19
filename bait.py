import HTSeq

gff_file = HTSeq.GFF_Reader( "/home/shashank/Downloads/gencode.v35.annotation.gtf.gz", end_included=True )

transcripts= {}

for feature in gff_file:
   if feature.type == "exon":
      transcript_id = feature.attr['gene_name']
      if transcript_id not in transcripts:
         transcripts[ transcript_id ] = list()
      transcripts[ transcript_id ].append( feature )
      
for transcript_id in sorted( transcripts ):      
   transcript_length = 0
   transcript_start = 0
   transcript_end = 0
   for exon in transcripts[ transcript_id ]:
      transcript_length += exon.iv.length + 1
      transcript_start = exon.iv.start 
      transcript_end = exon.iv.end
      chrsm = exon.iv.chrom
   print (chrsm, transcript_start,transcript_end,transcript_id, file=outF )
 
 outF = open("trans.bedpe","w")
