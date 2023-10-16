#!/usr/bin/env python3

"""
Author: Yiyun Li
Task: filter the network
"""

from sys import argv

def parse_network_file(filename):
    lines = open(filename).readlines()
    network = [] # network = [[regulator, target]]
    for line in lines:
        interaction = line.strip()
        interaction = interaction.split("\t")[0:2]
        network.append(interaction)
    return network

def get_regulator(gene, network):
    record = []
    for interaction in network:
        if gene == interaction[1]:
            record.append(interaction)
    return record

def get_2_regulator(gene, network):
    record = []
    for interaction in network:
        if gene == interaction[1]:
            record.append(interaction)
    re = [a[0] for a in record]
    for i in re:
        for interaction in network:
            if i == interaction[1]:
                if interaction not in record:
                    record.append(interaction)
    return record



def get_output(record, filename):
    with open(filename, 'w') as f:
        f.write("regulator\ttarget\n")
        for i in record:
            f.write("{}\t{}\n".format(i[0], i[1]))
    f.close()

def get_2_output(gene, record, filename):
    with open(filename, 'w') as f:
        f.write("regulator\ttarget\tlayer\n")
        for i in record:
            if gene in i[1]:
                i.append(1)
            else:
                i.append(2)

        for i in record:
            f.write("{}\t{}\t{}\n".format(i[0], i[1],i[2]))
    f.close()

def main():
    network_fn = argv[1]
    gene = argv[2]
    output_fn = argv[3]
    network_type = argv[4]
    network = parse_network_file(network_fn)
    if network_type == "1":
        record = get_regulator(gene, network)
        get_output(record, output_fn)
    if network_type == "2":
        record = get_2_regulator(gene, network)
        get_2_output(gene, record, output_fn)
    return

if __name__ == '__main__':
    main()