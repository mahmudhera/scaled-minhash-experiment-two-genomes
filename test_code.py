import subprocess
from compare_two_genomes import compare_two_files_to_get_multiple_containments

if __name__ == "__main__":
    f1 = "original.fasta"
    f2 = "mutated.fasta"
    num_runs = 100
    s = 0.1
    k = 21
    size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes = compare_two_files_to_get_multiple_containments(f1, f2, k, s, num_runs)
    
    print(scaled_containments)
    
    f = open("script.sh", 'w')
    for i in range(num_runs):
        seed_for_mash = i + 1
        command = "mash dist " + f1 + " " + f2 + " -s " +str(sketch_sizes[i])+ " -S " + str(seed_for_mash)
        f.write(command)
        f.write("\n")
    f.close()
    
    f = open('mash_output', 'w')
    cmd = "bash script.sh"
    cmd_args = cmd.split(' ')
    subprocess.call(cmd_args, stdout=f)
    f.close()
    
    f = open('mash_jaccards', 'w')
    cmd = 'cut -f3 mash_output'
    cmd_args = cmd.split(' ')
    subprocess.call(cmd_args, stdout=f)
    f.close()
    
    mash_jaccards = []
    f = open('mash_jaccards', 'r')
    lines = f.readlines()
    for line in lines:
        mash_jaccards.append( float(line.strip()) )
    print(mash_jaccards)
    f.close()
    
    mash_containments = []
    for j in mash_jaccards:
        c = j * 1.0 * size_union / size_1
        mash_containments.append(c)
    print(mash_containments)
    print(scaled_containments)
    