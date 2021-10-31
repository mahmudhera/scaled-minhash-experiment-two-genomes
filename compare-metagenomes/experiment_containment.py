import screed
import subprocess
import mmh3

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N"}

def reverse_complement(s):
    """
    Return reverse complement of 's'.
    """
    c = "".join(reversed([__complementTranslation[n] for n in s]))
    return c

def canonical_kmers(seq, k):
    for start in range(len(seq) - k + 1):
        kmer = seq[start: start + k].upper()
        rev_kmer = reverse_complement(kmer)
        if rev_kmer < kmer:
            kmer = rev_kmer

        if any(c not in "ACGT" for c in kmer):
            continue
        yield kmer

def get_kmers_in_file(filename, k):
	with screed.open(filename) as f:
		for record in f:
			for kmer in canonical_kmers(record.sequence, k):
				yield kmer

def get_hash_from_kmer(kmer, seed=0):
	hash_value = mmh3.hash64(kmer, seed=seed)[0]
	if hash_value < 0:
		hash_value += 2**64
	return hash_value

# c is float, 0 < c < 1
def extract_part_of_genome(c, genome_filename, out_filename):
    length = 0
    with screed.open(genome_filename) as f:
        for record in f:
            length += len(record.sequence)
    required_length = int(length * c)
    with screed.open(genome_filename) as f:
        for record in f:
            if len(record.sequence) < required_length:
                small_str = record.sequence[:required_length]
                break
    f = open(out_filename, 'w')
    f.write('> small_seq\n')
    f.write(small_str)
    f.close()
    
def create_super_metagenome(metagenome_filename, small_genome_filename, super_mg_filename):
    args = ['cat', metagenome_filename, small_genome_filename]
    f = open(super_mg_filename, 'w')
    subprocess.call(args, stdout=f)
    f.close()

def count_num_kmers_in_file(filename, k):
    kmer_hashes = set()
    for kmer in get_kmers_in_file(filename, k):
        kmer_hashes.add(get_hash_from_kmer(kmer))
    return len(kmer_hashes)

g_filename = 'ecoli.fasta'
smallg_filename = 'temp.fasta'
smg_filename = 'supermg.fasta'
mg_filename = 'SRR492065.contigs.fa'

scale_factor = 0.1
k = 21

extract_part_of_genome(0.1, g_filename, smallg_filename)
create_super_metagenome(mg_filename, smallg_filename, smg_filename)
num_kmers = count_num_kmers_in_file(g_filename, k)
expected_sketch_size = int(num_kmers * scale_factor)
print(expected_sketch_size, num_kmers)