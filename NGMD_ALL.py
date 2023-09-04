#!/opt/conda/envs/get/bin/python
import os, re, sys, subprocess, argparse, glob, shutil
from datetime import date
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_hybpiper_assemble(args):
    sample_names = []
    with open(args.hybpiper_config_file, 'r') as fi:
        for line in fi:
            infos1 = line.strip().split('\t')
            a = infos1[0]
            b = glob.glob(infos1[1])
            c = infos1[2]
            sample_names.append(c)

            for each_b in b:
                command = ["hybpiper", "assemble", "-t_dna", a, "-r", each_b, "--prefix", c, "--bwa", "--cpu", "10"]
                print(' '.join(command))
                subprocess.run(command)

    with open("sample_names_file.txt", "w") as out:
        for sample_name in sample_names:
            out.write(sample_name + "\n")

    args.sample_names_file = "sample_names_file.txt"

    with open(args.hybpiper_config_file, 'r') as fi:
        first_line = next(fi).strip().split('\t')
        args.target_file = first_line[0]

def run_hybpiper_retrieve(args):
    command = ["hybpiper", "retrieve_sequences", "-t_dna", args.target_file, "dna", "--sample_names", args.sample_names_file]
    print(' '.join(command))
    subprocess.run(command)

def prepare_for_seq_align_variant_call_hybpiper(args):
    fna_files = glob.glob(os.path.join(args.output_dir, '*.FNA'))
    return fna_files  # return list of .FNA files

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

def prepare_for_seq_align_variant_call_getorganelle(args):
    fasta_files = []
    all_sample_fasta = open("all_sample.fasta", "w")  # Output file to combine all fasta files
    with open(args.getorganelle_config_file, 'r') as fi:
        for line in fi:
            infos1 = line.strip().split('\t')
            output_dir = infos1[2]

            files = glob.glob(os.path.join(output_dir, '*.fasta'))
            for f in files:
                new_name = output_dir + "_" + os.path.basename(f)
                shutil.copy2(f, new_name)
                fasta_files.append(new_name)

                # Rename sequence ids in fasta file and write to the combined fasta file
                records = list(SeqIO.parse(f, "fasta"))
                for rec in records:
                    rec.id = output_dir
                    rec.description = output_dir
                    SeqIO.write(rec, all_sample_fasta, "fasta")
    
    all_sample_fasta.close()  # Always remember to close the file when you're done with it!
    return ["all_sample.fasta"]

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


def run_sequence_alignment_and_variant_calling(fasta_list):
    files_for_flanking_extraction = []
    for fasta_file in fasta_list:
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
		
        # Save the file names for flanking sequence extraction
        msafile = prefix+"_clustal.fasta"
        consensusfile = prefix+"_temp_cons.fasta"
        variantsfile = prefix+"_output.vcf"
        files_for_flanking_extraction.append((msafile, consensusfile, variantsfile))

    # Now, write the file names to a list file and run flanking sequence extraction
    list_file = "file_list.txt"
    with open(list_file, "w") as out_file:
        for msafile, consensusfile, variantsfile in files_for_flanking_extraction:
            out_file.write("\t".join([msafile, consensusfile, variantsfile]) + "\n")

    flanking_args = argparse.Namespace()
    flanking_args.files_list = list_file
    run_extract_flanking_sequences(flanking_args)
		
		

def run_extract_flanking_sequences(args):
    flanking_length = 300

    with open(args.files_list, 'r') as f:
        for line in f:
            files = line.strip().split('\t')
            msafile = files[0]
            consensusfile = files[1]
            variantsfile = files[2]

            prefix = os.path.splitext(msafile)[0].replace('_clustal', '')

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



def main():
    parser = argparse.ArgumentParser(description='Bioinformatics pipelines')
    subparsers = parser.add_subparsers(help='Choose one of the pipelines')

    parser_hybpiper = subparsers.add_parser('run_hybpiper_assemble_and_retrieve')
    parser_hybpiper.add_argument('hybpiper_config_file', help='config file for run_hybpiper_assemble')
    parser_hybpiper.add_argument('output_dir', help='Output directory for saving results')
    parser_hybpiper.set_defaults(func=lambda args: [run_hybpiper_assemble(args), run_hybpiper_retrieve(args), run_sequence_alignment_and_variant_calling(prepare_for_seq_align_variant_call_hybpiper(args))])

    parser_getorganelle = subparsers.add_parser('run_getorganelle')
    parser_getorganelle.add_argument('getorganelle_config_file', help='config file for run_getorganelle')
    parser_getorganelle.add_argument('target_database', help='target database for run_getorganelle:embplant_nr,embplant_pt,embplant_mt')
    parser_getorganelle.add_argument('output_dir', help='Output directory for saving results')
    parser_getorganelle.set_defaults(func=lambda args: [run_getorganelle(args), run_sequence_alignment_and_variant_calling(prepare_for_seq_align_variant_call_getorganelle(args))])

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

