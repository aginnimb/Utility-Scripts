import csv

'''
Description: This script reads the matrix generated from the hash-cgmlst tool and calculates the distances between two samples
date: 07/25/2023
Author: Anusha Ginni, qxu0@cdc.gov
usage:  python hash-cgmlst_dist2.py. Dont forget change the path to the tsv file with in the script'''

def read_tsvfile(file_path):
    data = {}
    with open(file_path, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        header = next(reader)
        for row in reader:
            sample_id = row[0] # assembly Ids
            alleles = [(allele) for allele in row[1:]] #hashvalues of allele numbers
            # print(alleles)
            data[sample_id] = alleles
            # print(data)
    return data

def calculate_distance(sample1, sample2):
    total_loci = len(sample1)
    common_alleles = sum(1 for allele1, allele2 in zip(sample1, sample2) if allele1 == allele2 != "-1") #get common alleles between two samples
    # same_loci = sum(1 for allele1, allele2 in zip(sample1, sample2) if allele1 == allele2 != 0)
    #excluding loci with -1 from the compared
    allele_compared = sum(1 for allele1, allele2 in zip(sample1, sample2) if allele1 != "-1" or allele2 != "-1")
    percent_identity = (common_alleles / allele_compared) * 100
    return percent_identity, common_alleles, allele_compared, total_loci

# def calculate_distance(sample1, sample2):
#     total_loci = len(sample1)
#     common_alleles = 0
#     same_loci = 0
#     allele_compared = 0
#     exclude_allele1_1 = sample1.count("-1")
#     exclude_allele_2 = sample2.count("-1")
#     allele_excluded = int(exclude_allele1_1) + int(exclude_allele_2)
#     for allele1, allele2 in zip(sample1, sample2):
#         if allele1 == -1 or allele2 == -1:
#             continue
#         elif allele1 == allele2 != -1:
#             common_alleles += 1 #get common alleles between two samples
#         elif allele1 != 1 or allele2 != 1:
#             allele_compared += 1
#         else:
#             same_loci += 1
#     percent_identity = (common_alleles / allele_compared) * 100
#     return percent_identity, common_alleles, allele_compared, total_loci

def write_csv_output(output_file, data):
    with open(output_file, 'w', newline='') as csvfile:
        colnames = ['sample1', 'sample2', 'identity', 'num_alleles_same','num_allele_compared','total_loci']
        writer = csv.DictWriter(csvfile, fieldnames=colnames)
        writer.writeheader()
        samples = list(data.keys())
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                sample1 = samples[i]
                sample2 = samples[j]
                percent_identity,num_alleles_same,num_allele_compared,total_loci = calculate_distance(data[sample1], data[sample2])
                writer.writerow({'sample1': sample1, 'sample2': sample2, 'identity': percent_identity,'num_alleles_same': num_alleles_same, 'num_allele_compared': num_allele_compared, 'total_loci': total_loci})

if __name__ == '__main__':
    input_tsv_file = 'hash-cgmlst.matrixupdated.tsv'
    output_csv_file = 'hash-cgmlst.matrixupdated.csv'
    data = read_tsvfile(input_tsv_file)
    # print(data)
    write_csv_output(output_csv_file, data)