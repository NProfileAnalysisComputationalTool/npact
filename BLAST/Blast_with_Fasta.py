#blast: blastp
#database:  nr
#usage: python Blast_with_Fasta.py fasta_file e_threshold > output_file
#output:SEQUENCE ID, SIZE, NUMBER OF HITS: 
#       alignment title(description), alignemnt length, e-value, score, identity, accession

import sys

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

fasta_file = sys.argv[1]
E_VALUE_THRESH = float(sys.argv[2])

output_file = fasta_file+".blast"
for seq_record in SeqIO.parse(fasta_file, "fasta"):
	print "SEQUENCE ID:",seq_record.id,
	print "SIZE:", len(seq_record)
	fasta_string= seq_record.seq
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)#,format_type="Text")
	blast_records = NCBIXML.parse(result_handle)
	for blast_record in blast_records:
		print "NUMBER OF HITS:", len(blast_record.alignments)
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					print alignment.hit_def+"\t"+str(alignment.length)+"\t"+str(hsp.expect)+"\t"+str(hsp.score)+"\t"+str(hsp.identities)+"\t"+str(alignment.accession)




# for more informations..
 				print 'e-value:', hsp.expect
 				print "score", hsp.score
 				print "accession", alignment.accession
				print "query_start_end",hsp.query_start, hsp.query_end
				print "subject_start_end",hsp.sbjct_start, hsp.sbjct_end
				print "identity", hsp.identities
				print "query:",hsp.query
 				print "match:",hsp.match
 				print "target:",hsp.sbjct
 				
