from compare_two_genomes import compare_two_files_to_get_multiple_containments

if __name__ == "__main__":
    f1 = "original.fasta"
    f2 = "mutated.fasta"
    num_runs = 100
    s = 0.1
    k = 21
    values = compare_two_files_to_get_multiple_containments(f1, f2, k, s, num_runs)
    print(values)