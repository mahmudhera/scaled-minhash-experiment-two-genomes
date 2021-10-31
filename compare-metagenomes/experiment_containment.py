import screed
import subprocess
import mmh3

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N"}

class ScaledMinHash:
    def __init__(self, scale_factor, max_hash_value):
        self.hash_set = set()
        self.H = max_hash_value
        self.scale_factor = scale_factor
        self.raw_elements = set()
        
    def add_value(self, hash_value):
        if hash_value <= self.H * self.scale_factor:
            self.hash_set.add(hash_value)
        self.raw_elements.add(hash_value)
            
    def add_values(self, hash_values):
        for hash_value in hash_values:
            self.add_value(hash_value)
            
    def remove(self, hash_value):
        self.hash_set -= hash_value
        
    def print_hash_set(self):
        print(self.H, self.scale_factor, self.hash_set)
        
    def get_containment(self, smh):
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / len(self.hash_set)
    
    def get_scaled_containment(self, smh):
        bf = 1 - (1 - self.scale_factor) ** len(self.raw_elements)
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / ( len(self.hash_set) * bf )

    def get_sketch_size(self):
        return len( self.hash_set )

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

def compare_two_files_to_get_multiple_containments(filename_1, filename_2, k, scale_facor, num_runs):
    seeds = [i+1 for i in range(num_runs)]
    H = int(2**64)

    sketch_sizes = []
    scaled_containments = []
    for seed in seeds:
        kmer_hashes_1 = set()
        kmer_hashes_2 = set()
        for kmer in get_kmers_in_file(filename_1, k):
            kmer_hashes_1.add(get_hash_from_kmer(kmer, seed=seed))
        for kmer in get_kmers_in_file(filename_2, k):
            kmer_hashes_2.add(get_hash_from_kmer(kmer, seed=seed))
            
        size_1 = len(kmer_hashes_1)
        size_2 = len(kmer_hashes_2)
        size_union = len( kmer_hashes_1.union( kmer_hashes_2 ) )
        size_intersection = len( kmer_hashes_1.intersection( kmer_hashes_2 ) )            
        
        smh1 = ScaledMinHash(scale_facor, H)
        smh1.add_values(kmer_hashes_1)
        smh2 = ScaledMinHash(scale_facor, H)
        smh2.add_values(kmer_hashes_2)
        
        scaled_containment = smh1.get_containment(smh2)
        sketch_size = smh1.get_sketch_size()
        
        sketch_sizes.append(sketch_size)
        scaled_containments.append(scaled_containment)
    
    true_containment = 1.0*size_intersection/size_1
    return size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes
    

g_filename = 'ecoli.fasta'
smallg_filename = 'temp.fasta'
smg_filename = 'supermg.fasta'
mg_filename = 'SRR492065.contigs.fa'

scale_factor = 0.002
k = 21
C = 0.1
num_runs = 10
seed = 1

extract_part_of_genome(C, g_filename, smallg_filename)
create_super_metagenome(mg_filename, smallg_filename, smg_filename)
num_kmers = count_num_kmers_in_file(g_filename, k)
expected_sketch_size = int(num_kmers * scale_factor)
size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes = compare_two_files_to_get_multiple_containments(g_filename, smg_filename, k, scale_factor, num_runs)
print(true_containment)
print(scaled_containments)