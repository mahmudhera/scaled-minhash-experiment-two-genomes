import screed
import subprocess

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


g_filename = 'ecoli.fasta'
smallg_filename = 'temp.fasta'
smg_filename = 'supermg.fasta'
mg_filename = 'SRR492065.contigs.fa'

extract_part_of_genome(0.1, g_filename, smallg_filename)
create_super_metagenome(mg_filename, smallg_filename, smg_filename)
