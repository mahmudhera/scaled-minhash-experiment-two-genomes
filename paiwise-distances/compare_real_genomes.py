import screed
import glob
import mmh3
import subprocess
import string

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N", 'K': 'N'}
for char in string.ascii_uppercase:
    if char not in __complementTranslation.keys():
        __complementTranslation[char] = 'N'

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

def count_num_kmers_in_file(filename, k):
    kmer_hashes = set()
    for kmer in get_kmers_in_file(filename, k):
        kmer_hashes.add(get_hash_from_kmer(kmer))
    return len(kmer_hashes)

def get_true_mut_rate(filename1, filename2):
    cmd = "java -jar OAU.jar -u ./usearch11.0.667_i86linux32 --f1 " + filename1 + " -f2 " + filename2 + " ."
    args = cmd.split(' ')
    f = open('temp', 'w')
    subprocess.call(args, stdout=f)
    f.close()
    f = open('temp', 'r')
    true_ani = float(f.readlines()[-1].split('\t')[1])
    f.close()
    return 100.0-true_ani
    

seed = 2
stats_filename = 'results'
k = 21
scale_factor = 0.01
num_runs = 2

f = open(stats_filename, 'w')
f.close()

for filename1 in glob.glob('*.fna'):
    for filename2 in glob.glob('*.fna'):
        if filename1 == filename2:
            continue
    
    mutation_rate = get_true_mut_rate(filename1, filename2)
    cmd = "python test_code.py " + filename1 + " " + filename2 + " -k " + str(k) + " -s " + str(scale_factor) + " --seed " + str(seed) + " -c 0.95 -N " + str(num_runs) + " -p " + str(mutation_rate) + " --fout " + stats_filename
    args = cmd.split(' ')
    subprocess.call(args)
    
