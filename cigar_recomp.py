import pysam
import sys

infile = sys.argv[1]
fa = sys.argv[2]
outfile = sys.argv[3]

bamfile = pysam.AlignmentFile(infile, 'rb')
fafile = pysam.FastaFile(fa)

fout = pysam.AlignmentFile(outfile, 'w', header=bamfile.header)

def rc_seq(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return "".join(complement.get(base, base) for base in reversed(seq))

for read in bamfile.fetch():
	if read.cigarstring is not None:
			cigar = read.cigartuples
			ref_start = read.reference_start
			ref_end = read.reference_end
			ref_id = read.reference_name
			qry_start = read.query_alignment_start
			qry_end = read.query_alignment_end
			qry_aln_len = read.query_alignment_length
			qry_name = read.query_name
			newCigar = ""
			rOff, qOff = 0, 0
			if read.is_reverse:
				qry_seq = rc_seq(fafile.fetch(qry_name))
			else:
				qry_seq = fafile.fetch(qry_name, qry_start, qry_end)
			ref_seq = fafile.fetch(ref_id, ref_start, ref_end)
			#print "ref_start", ref_start, "qry_start", qry_start, "reverse", read.is_reverse
			#print ref_seq
			#print
			#print qry_seq
			nm = 0
			for (operation, operation_len) in cigar:
					#print "operation", operation, "operation_len", operation_len
					#print "rOff", rOff, "qOff", qOff
					if operation == 0:
							#print "qry match part", qry_seq[qOff:qOff+operation_len], len(qry_seq[qOff:qOff+operation_len])
							#print "ref match part", ref_seq[rOff:rOff+operation_len], len(ref_seq[rOff:rOff+operation_len])
							pre_state = None
							state_len = 0
							for i in range(operation_len):
									if ref_seq[rOff+i] == qry_seq[qOff+i]:
											state = '='
									else:
											state = 'X'
									if pre_state is None:
											pre_state = state
											state_len += 1
											if state == 'X':
												nm += 1
									else:
											if pre_state != state:
													newCigar += "%d%s" %(state_len, pre_state)
													pre_state = state
													state_len = 1
													if state == 'X':
														nm += 1
											else:
													state_len += 1
													if state == 'X':
														nm += 1
							newCigar += "%d%s" %(state_len, pre_state)
							rOff += operation_len
							qOff += operation_len
							#print "M Intercept", "newCigar", newCigar
							#print "NM", nm
					else:
							if operation == 4:
									newCigar += "%dS" %(operation_len)
							if operation == 1:
									newCigar += "%dI" %(operation_len)
									qOff += operation_len
							if operation == 2:
									newCigar += "%dD" %(operation_len)
									rOff += operation_len
							if operation == 5:
									newCigar += "%dH" %(operation_len)
			#print "finish NM", nm
			read.cigarstring = newCigar
			fout.write(read)
	else:
		# just write out
		fout.write(read)
fout.close()
fafile.close()
bamfile.close()


