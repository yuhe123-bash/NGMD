#!/opt/conda/envs/get/bin/python
import os, re, sys, subprocess, argparse, glob
from datetime import date
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_hybpiper_assemble(args):
    with open(args.hybpiper_assemble_config_file, 'r') as fi:
        for line in fi:
            infos1 = line.strip().split('\t')
            a = infos1[0]
            b = glob.glob(infos1[1])
            c = infos1[2]

            for each_b in b:
                command = ["hybpiper", "assemble", "-t_dna", a, "-r", each_b, "--prefix", c, "--bwa", "--cpu", "10"]
                print(' '.join(command))
                subprocess.run(command)

def run_hybpiper_retrieve(args):
    command = ["hybpiper", "retrieve_sequences", "-t_dna", args.target_file, "dna", "--sample_names", args.sample_names_file]
    print(' '.join(command))
    subprocess.run(command)

def run_getorganelle(args):
    with open(args.getorganelle_config_file, 'r') as fi:
        for line in fi:
            infos1 = line.strip().split('\t')
            a = infos1[0]
            b = infos1[1]
            c = infos1[2]

            if args.target_database == "embplant_nr":
                command = ["get_organelle_from_reads.py", "-1", a, "-2", b, "-o", c, "-R", "10", "-k", "35,85,115", "-F", args.target_database]
            elif args.target_database == "embplant_pt":
                command = ["get_organelle_from_reads.py", "-1", a, "-2", b, "-o", c, "-R", "15", "-k", "21,45,65,85,105", "-F", args.target_database]
            elif args.target_database == "embplant_mt":
                command = ["get_organelle_from_reads.py", "-1", a, "-2", b, "-o", c, "-R", "20", "-k", "21,45,65,85,105", "-P", "1000000", "-F", args.target_database]
            else:
                print(f"Invalid target_database: {args.target_database}. Skipping.")
                continue

            print(' '.join(command))
            subprocess.run(command)

def run_clustalw2(input_file, prefix):
    clustalw_exe = "~/.conda/envs/hybpiper/bin/clustalw2"  # Path to ClustalW2 executable
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=input_file, outfile=prefix+"_clustal.fasta", output="fasta")
    stdout, stderr = clustalw_cline()

    align = AlignIO.read(prefix+"_clustal.fasta", "fasta")
    return align

def get_consensus_sequence(align, prefix):
    tmp_file = prefix+"_temp.fasta"
    AlignIO.write(align, tmp_file, "fasta")
    consensus_command = "cons -sequence " + tmp_file + " -outseq " + prefix + "_temp_cons.fasta"
    process = subprocess.Popen(consensus_command, shell=True)
    process.wait()
    os.remove(tmp_file)
    with open(prefix+"_temp_cons.fasta") as cons_file:
        for record in SeqIO.parse(cons_file, "fasta"):
            consensus_sequence = str(record.seq)
    return consensus_sequence

def run_sequence_alignment_and_variant_calling(args):
    with open(args.fasta_list, 'r') as fi:
        for line in fi:
            fasta_file = line.strip()
            prefix = os.path.splitext(fasta_file)[0]
            align = run_clustalw2(fasta_file, prefix)
            refseq = get_consensus_sequence(align, prefix)

            # Write the snptable.txt file
            snp_output_file = open(prefix + "_snptable.txt", "w")
            for record in align:
                snp_output_file.write(record.id + "\t" + str(record.seq) + "\n")
            snp_output_file.close()

            version = "1.0"  # Replace with your version
            # vcf header
            vcfheader = '''##fileformat=VCFv4.2
            ##contig=<ID=mychrom,length=chromlength>
            ##fileDate=todaydate
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##source=getvar.py
            ##getvar.pyVersion="myversion"
            ##getvar.pyCmd="mycmd"
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'''
            vcfheader = vcfheader.replace("myversion", version)
            # add date
            today = date.today()
            d1 = today.strftime("%Y%m%d")
            vcfheader = vcfheader.replace("todaydate", d1)
            refseq_nogap = refseq.replace("-","")
            refseq_length = len(refseq_nogap)
            # update header
            vcfheader = vcfheader.replace("mychrom", "consensus")
            vcfheader = vcfheader.replace("chromlength", str(refseq_length))
            seqnames = [rec.id for rec in align]
            vcfheader = vcfheader + "\t" + "\t".join(seqnames)

            vcf_output_file = open(prefix + "_output.vcf", "w")
            vcf_output_file.write(vcfheader + "\n")

            n_nogap = -1 # no gap position for reference
            n = 0 # no gap position for mutations
            seq2 = {} # new dictionary of rotated sequences with no-gap postions as keys
            for i in range(len(refseq)): # i is the position of seq position
                if refseq[i] != "-":
                    n_nogap += 1
                templist = []
                nblank = 0
                for seq_record in align:
                    seq = str(seq_record.seq)
                    templist.append(seq[i])
                    if seq[i] == "-":
                        nblank += 1
                if nblank == 0:
                    n = n_nogap
                if n not in seq2:
                    seq2[n] = templist
                else:
                    seq2[n] = [s1 + s2 for s1, s2 in zip(seq2[n], templist)]

            outlist = []
            for k in range(len(refseq_nogap)):
                if (k in seq2):
                    seq = [w.replace('-', '') if len(w) > 1 else w for w in seq2[k]]
                    alleles = set(seq)
                    if len(alleles) != 1:
                        ref_allele = seq[0]
                        alt_allele_set = alleles - set([ref_allele])
                        alt_allele = ",".join(alt_allele_set)
                        outlist.append([str(k+1), ref_allele, alt_allele] + seq)

            # write vcf to a file
            for i in outlist:
                alist = i[1:2] + i[2].split(",") # allele list
                indexlist = [str(alist.index(x))+"/" + str(alist.index(x)) for x in i[3:]]
                j = ["consensus"] + i[0:3] + [".", ".", ".", "GT"] + indexlist
                vcf_output_file.write("\t".join(j) + "\n")
            vcf_output_file.close()


def run_extract_flanking_sequences(args):
    flanking_length = 300

    with open(args.files_list, 'r') as f:
        for line in f:
            files = line.strip().split('\t')
            msafile = files[0]
            consensusfile = files[1]
            variantsfile = files[2]
            
            # Extract the prefix from msafile to be used in the output filename
            prefix = os.path.splitext(msafile)[0]

            fasta_records = SeqIO.to_dict(SeqIO.parse(msafile, "fasta"))
            consensus_record = SeqIO.read(consensusfile, "fasta")

            with open(variantsfile) as vf:
                for line in vf:
                    if line.lstrip().startswith("#CHROM"):
                        individual_ids = line.lstrip().split("\t")[9:]
                        break

            try:
                individual_ids
            except NameError:
                print("Error: Could not parse individual IDs from VCF file header.")
                sys.exit(1)

            new_records = []

            with open(variantsfile) as vf:
                for line in vf:
                    if not line.startswith("#"):
                        items = line.strip().split("\t")
                        if len(items) < 2 or items[0] == "#CHROM":
                            continue
                        variant_pos = int(items[1])
                        genotypes = items[9:]

                        pos_in_msa = {pos: i for pos, i in enumerate(j for j, base in enumerate(str(consensus_record.seq)) if base != '-')}

                        for genotype, individual_id in zip(genotypes, individual_ids):
                            if genotype not in ["0/0", ".", "./.", "0|0", ".|."]:
                                record = fasta_records[individual_id]
                                msa_variant_pos = pos_in_msa[variant_pos-1]
                                if msa_variant_pos - flanking_length >= 0 and msa_variant_pos + flanking_length < len(record.seq):
                                    flanking_seq = str(record.seq[msa_variant_pos-flanking_length:msa_variant_pos+flanking_length])
                                    flanking_seq = flanking_seq.replace("-", "")
                                    new_record = SeqRecord(Seq(flanking_seq), id=f"{individual_id}_pos_{variant_pos}", description="")
                                    new_records.append(new_record)

            # Generate a unique filename for each run
            output_filename = f"{prefix}_flanking_sequences.fasta"
            SeqIO.write(new_records, output_filename, "fasta")





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    parser_hybpiper = subparsers.add_parser('run_hybpiper_assemble')
    parser_hybpiper.add_argument('hybpiper_assemble_config_file', help='config file for run_hybpiper_assemble')
    parser_hybpiper.set_defaults(func=run_hybpiper_assemble)

    parser_hybpiper_retrieve = subparsers.add_parser('run_hybpiper_retrieve')
    parser_hybpiper_retrieve.add_argument('sample_names_file', help='sample names file for run_hybpiper_retrieve')
    parser_hybpiper_retrieve.add_argument('target_file', help='target file for run_hybpiper_retrieve')
    parser_hybpiper_retrieve.set_defaults(func=run_hybpiper_retrieve)

    parser_getorganelle = subparsers.add_parser('run_getorganelle')
    parser_getorganelle.add_argument('getorganelle_config_file', help='config file for run_getorganelle')
    parser_getorganelle.add_argument('target_database', help='target database for run_getorganelle')
    parser_getorganelle.set_defaults(func=run_getorganelle)
	
    parser_variant_calling = subparsers.add_parser('run_variant_calling')
    parser_variant_calling.add_argument('fasta_list', help='file containing list of fasta files for variant calling')
    parser_variant_calling.set_defaults(func=run_sequence_alignment_and_variant_calling)


    parser_flanking = subparsers.add_parser('run_extract_flanking_sequences')
    parser_flanking.add_argument('files_list', help='file containing msafile, consensusfile and variantsfile for each individual')
    parser_flanking.set_defaults(func=run_extract_flanking_sequences)
	
    args = parser.parse_args()
    args.func(args)