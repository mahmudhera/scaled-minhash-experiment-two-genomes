import screed

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
    
extract_part_of_genome(0.1, 'ecoli.fasta', 'test.fasta')