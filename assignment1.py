#!/usr/bin/python3

import mysql.connector
import pysam
import os

__author__ = 'Richard Kriz'


class Assignment1:

    def __init__(self, gene, bam_file, coord_file):
        self.gene = gene
        self.bam_file = bam_file
        self.coord_file = coord_file
        self.sam_file = pysam.AlignmentFile(self.bam_file, "rb")

    def download_gene_coordinates(self, genome_reference):
        print("Connecting to UCSC to fetch data")

        # Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep',
                                      passwd='password', db=genome_reference)

        # Get cursor
        cursor = cnx.cursor()

        # Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        # Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

        # Execute query
        cursor.execute(query)
        with open(self.coord_file, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")

        # Close cursor & connection
        cursor.close()
        cnx.close()

        print("Done fetching data")

    def get_gene_specifications(self):
        with open(self.coord_file) as coord_file:
            for row in coord_file:
                if self.gene in row:
                    row = row[1:-2].replace("'", "").replace(" ", "")
                    row = row.split(",")
                    break
        self.chr = row[2]
        self.start = int(row[3])
        self.stop = int(row[4])
        self.exons = int(row[6])
        exon_regions_delimiter = int((len(row)-7)/2+6)      # position of delimiter for start and stop coords for exons
        exon_regions = row[7:exon_regions_delimiter] + row[exon_regions_delimiter+1:-1]
        self.exon_coord = []
        for i in range(int(len(exon_regions)/2)):
            self.exon_coord.append(exon_regions[i] + " - " + exon_regions[int(len(exon_regions)/2+i)])

    def get_sam_header(self):
        self.sam_header = {}
        for key, value in self.sam_file.header["HD"].items():
            if key == "VN":
                self.sam_header[0] = value
            if key == "SO":
                self.sam_header[1] = value

    def get_properly_paired_reads_of_gene(self):
        self.reads = list(self.sam_file.fetch(self.chr, self.start, self.stop))
        self.reads_proper =len([i for i in self.reads if i.is_proper_pair])

        return self.reads_proper

    def get_gene_reads_with_indels(self):
        self.reads_indel = 0
        for i in self.reads:
            if not i.is_unmapped:
                cig = i.cigartuples
                for operator, length in cig:
                    if (operator == 1) or (operator == 2):
                        self.reads_indel += 1
        return self.reads_indel

    def calculate_total_average_coverage(self):
        for i in self.sam_file.header["SQ"]:
            if i["SN"] == self.chr:
                self.chr_len = i["LN"]

        # ToDo: comment out if testing
        coverage_matrix = self.sam_file.count_coverage(self.chr, start=0, stop=self.chr_len)
        coverage_list_sum = 0
        for i in range(len(coverage_matrix)):
            for j in coverage_matrix[0]:
                coverage_list_sum += coverage_matrix[i][j]
        self.total_cov = round(coverage_list_sum / len(coverage_matrix[0]), 2)

        return self.total_cov

    def calculate_gene_average_coverage(self):
        coverage_matrix = self.sam_file.count_coverage(self.chr, start=self.start, stop=self.stop)
        coverage_list_sum = 0
        for i in range(len(coverage_matrix)):
            for j in coverage_matrix[0]:
                coverage_list_sum += coverage_matrix[i][j]
        self.gene_cov = round(coverage_list_sum / len(coverage_matrix[0]), 2)

        return self.gene_cov

    def get_number_mapped_reads(self):
        self.reads_mapped = 0
        for i in self.reads:
            if not i.is_unmapped:
                self.reads_mapped += 1
        return self.reads_mapped

    def print_summary(self):
        self.get_gene_specifications()
        print("\nGene Symbol: {}".format(self.gene))
        print("Gene Coordinates: {} - {}".format(self.start, self.stop))
        self.get_sam_header()
        print("Sam Header:    VN: {}    Sorting: {}".format(self.sam_header[0], self.sam_header[1]))
        print("Number of properly paired reads: {}".format(self.get_properly_paired_reads_of_gene()))
        print("Gene Reads with INDELs: {}".format(self.get_gene_reads_with_indels()))
        print("Calculating Total Average Coverage. This can take a while...")
        print("Total Average Coverage: {}".format(self.calculate_total_average_coverage()))
        print("Gene Average Coverage: {}".format(self.calculate_gene_average_coverage()))
        print("Number of mapped reads: {}".format(self.get_number_mapped_reads()))
        print("Number of exons: {}".format(self.exons))
        print("Exonic Regions of Gene:")
        for i in self.exon_coord:
            print(i)


def main():
    print("Assignment 1")
    assignment1 = Assignment1("\'CBS\'", "./../chr21.bam", "./../hg19_data")
    if not os.path.isfile(assignment1.coord_file):
        assignment1.download_gene_coordinates("hg19")
    assignment1.print_summary()

    print("Done with assignment 1")


if __name__ == '__main__':
    main()
