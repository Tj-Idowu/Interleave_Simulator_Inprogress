#!/usr/bin/env python3

import os
import gzip

# Input path
genomes_folder = "/mnt/datafive/tj/datasets/Simulated_genomes"
genomes = os.listdir(genomes_folder)

# Output paths
bacillus = "/mnt/datafive/tj/datasets/Simulated_metagenomes/Simulated_WastewaterwoBacillus_R2.fastq.gz"
clostridium = "/mnt/datafive/tj/datasets/Simulated_metagenomes/Simulated_WastewaterwoClostridium_R2.fastq.gz"
escherichia = "/mnt/datafive/tj/datasets/Simulated_metagenomes/Simulated_WastewaterwoEscherichia_R2.fastq.gz"
francisella = "/mnt/datafive/tj/datasets/Simulated_metagenomes/Simulated_WastewaterwoFrancisella_R2.fastq.gz"
yersinia = "/mnt/datafive/tj/datasets/Simulated_metagenomes/Simulated_WastewaterwoYersinia_R2.fastq.gz"

# Make the Bacillus background with os module
with gzip.open(bacillus, 'wb') as outfile:
    for genome in os.listdir(genomes_folder):
        if genome.endswith("_R2.fastq.gz"):
            if genome.startswith("sbd_Bacillus_"):
                continue
            full_path = os.path.join(genomes_folder, genome)
            with gzip.open(full_path, 'rb') as infile:
                content = infile.read()
                outfile.write(content)
                if not content.endswith(b'\n'):
                    outfile.write(b'\n')

# Make the Clostridium background with os module
with gzip.open(clostridium, 'wb') as outfile:
    for genome in os.listdir(genomes_folder):
        if genome.endswith("_R2.fastq.gz"):
            if genome.startswith("sbd_Clostridium_"):
                continue
            full_path = os.path.join(genomes_folder, genome)
            with gzip.open(full_path, 'rb') as infile:
                content = infile.read()
                outfile.write(content)
                if not content.endswith(b'\n'):
                    outfile.write(b'\n')

# Make the Escherichia background with os module
with gzip.open(escherichia, 'wb') as outfile:
    for genome in os.listdir(genomes_folder):
        if genome.endswith("_R2.fastq.gz"):
            if genome.startswith("sbd_Escherichia_"):
                continue
            full_path = os.path.join(genomes_folder, genome)
            with gzip.open(full_path, 'rb') as infile:
                content = infile.read()
                outfile.write(content)
                if not content.endswith(b'\n'):
                    outfile.write(b'\n')

# Make the Francisella background with os module
with gzip.open(francisella, 'wb') as outfile:
    for genome in os.listdir(genomes_folder):
        if genome.endswith("_R2.fastq.gz"):
            if genome.startswith("sbd_Francisella_"):
                continue
            full_path = os.path.join(genomes_folder, genome)
            with gzip.open(full_path, 'rb') as infile:
                content = infile.read()
                outfile.write(content)
                if not content.endswith(b'\n'):
                    outfile.write(b'\n')

# Make the Yersinia background with os module
with gzip.open(yersinia, 'wb') as outfile:
    for genome in os.listdir(genomes_folder):
        if genome.endswith("_R2.fastq.gz"):
            if genome.startswith("sbd_Yersinia_"):
                continue
            full_path = os.path.join(genomes_folder, genome)
            with gzip.open(full_path, 'rb') as infile:
                content = infile.read()
                outfile.write(content)
                if not content.endswith(b'\n'):
                    outfile.write(b'\n')

