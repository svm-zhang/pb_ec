import pysam
import sys
import multiprocessing as mp
import os

def rc_seq(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return "".join(complement.get(base, base) for base in reversed(seq))

def worker(task_q, result_q, infile, fa, outp):
	while True:
		try:
			sub_refs, nth_job = task_q.get()
			refs, ref_lens = [], []
			for ref in sub_refs:
				refs += [ref[0]]
				ref_lens += [ref[1]]
			fafile = pysam.Fastafile(fa)
			outfile = outp + "%d.bam" %(nth_job)
			print "# refs", len(refs), "nth_job", nth_job, "outfile", outfile
			bamfile = pysam.AlignmentFile(infile, 'rb')
			#fout = pysam.AlignmentFile(outfile, 'w', reference_names=refs, reference_lengths=ref_lens)
			fout = pysam.AlignmentFile(outfile, 'wb', header=bamfile.header)
			for ref in refs:
				for read in bamfile.fetch(reference=ref, until_eof=True):
					if read.cigarstring is not None:
						cigar = read.cigartuples
						ref_start = read.reference_start
						ref_end = read.reference_end
						ref_id = read.reference_name
						qry_start = read.query_alignment_start
						qry_end = read.query_alignment_end
						qry_aln_len = read.query_alignment_length
						qry_name = read.query_name
						#newCigar = ""
						newCigar = ()
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
									else:
										if pre_state != state:
											#newCigar += "%d%s" %(state_len, pre_state)
											if pre_state == '=':
												newCigar += ((7, state_len), )
											elif pre_state == 'X':
												newCigar += ((8, state_len), )
											pre_state = state
											state_len = 1
										else:
											state_len += 1
								#newCigar += "%d%s" %(state_len, pre_state)
								if pre_state == '=':
									newCigar += ((7, state_len), )
								elif pre_state == 'X':
									newCigar += ((8, state_len), )
								rOff += operation_len
								qOff += operation_len
								#print "M Intercept", "newCigar", newCigar
							else:
								newCigar += ((operation, operation_len), )
						#read.cigarstring = newCigar
						read.cigar = newCigar
						fout.write(read)
					else:
						# just write out
						fout.write(read)
			bamfile.close()
			fafile.close()
			fout.close()
			result_q.put(outfile)
			print nth_job, "DONE"
		finally:
			task_q.task_done()

def main():
	infile = sys.argv[1]
	fa = sys.argv[2]
	outp = sys.argv[3]
	nproc = int(sys.argv[4])

	fafile = pysam.Fastafile(fa)
	nseq = len(fafile.references)
	nseq_per_proc = nseq/nproc + 1
	print "nseq", nseq, "nseq_per_proc", nseq_per_proc
	refs = []
	for i, length in enumerate(fafile.lengths):
		refs += [(fafile.references[i], length)]
	fafile.close()

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	for _ in range(nproc):
		p = mp.Process(target=worker, args=(task_q, result_q, infile, fa, outp))
		p.daemon = True
		p.start()

	i = 0
	nth_job = 1
	while i+nseq_per_proc <= nseq:
		task_q.put( (refs[i:i+nseq_per_proc+1], nth_job) )
		i += nseq_per_proc + 1
		nth_job += 1
	task_q.put( (refs[i:], nth_job) )

	try:
		task_q.join()
	except KeyboardInterrupt:
		sys.exit()
	else:
		bams = []
		while nproc:
			file = result_q.get()
			bams += [os.path.join(os.getcwd(), file)]
			print os.path.join(os.getcwd(), file)
			nproc -= 1
		#print bams
		#pysam.cat("-o", "%s.all.bam" %(outp), " ".join(bams))

if __name__ == "__main__":
	main()

