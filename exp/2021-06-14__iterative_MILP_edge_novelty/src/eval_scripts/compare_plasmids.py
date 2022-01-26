#!/usr/bin/python

# Analyses how well predicted plasmids resp. their contigs cover reference plasmids.
# It determines the proportions of reference plasmids covered by individual / all predicted plasmids (and vice versa)
# and uses these information to score the predictions (recall, precision, F1 score).
from __future__ import division
import argparse
import pandas as pd
from Bio import SeqIO



# read BLAST output file into table
def read_blast_output(file):
    col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]  # outfmt 6
    return pd.read_csv(file, sep = '\t', names = col_names, dtype = str)


# compute number of positions (integers) covered by a list of potentially overlapping intervals
def num_covered_positions(intervals):
    intervals.sort(key = lambda x: x[0])  # intervals is now sorted by start position

    num_pos_covered = 0
    last_pos_covered = 0  # last (right-most) position of contig covered so far
    for start, end in intervals:
        if end <= last_pos_covered:
            pass  # contained in previous interval -> no new position covered
        else:
            num_pos_covered += end - max(last_pos_covered + 1, start) + 1
            last_pos_covered = end

    return num_pos_covered


# convert information into a single row of a CSV-formatted file
def print_csv_row(tool, sample_id, values, sep):
    print(tool + sep + sample_id + sep + sep.join(map(str, values)))


# analyse matches between the predicted plasmids (= concatenation of orientated contigs) obtained from greedy strategy
# and the reference plasmids
def analyse_greedy(hits, predicted_sequences, reference_sequences, ref_chromosomes, ref_plasmids, threshold, name, sample_id, out):

    predicted_plasmids = sorted(list(predicted_sequences.keys()))
    covered_ref_sections = dict([(ref, []) for ref in reference_sequences.keys()]) # of lists
    covered_pred_sections = dict([(pred, []) for pred in predicted_plasmids])
    covered_ref_sections_per_pred = dict()
    covered_pred_sections_per_ref = dict()

    for pred in predicted_plasmids:
        for ref in sorted(ref_plasmids):
        #for ref in sorted(ref_chromosomes):
            pred_ref_hits = hits.loc[hits.qseqid == pred].loc[hits.sseqid == ref]
            covered_pred_sections_per_ref[(pred, ref)] = []
            covered_ref_sections_per_pred[(ref, pred)] = []

            for index, row in pred_ref_hits.iterrows():
                sstart = int(row[8])
                send = int(row[9])
                interval = (sstart, send) if sstart <= send else (send, sstart)
                covered_ref_sections[ref].append(interval)
                covered_ref_sections_per_pred[(ref, pred)].append(interval)

                qstart = int(row[6])
                qend = int(row[7])
                interval = (qstart, qend) if qstart <= qend else (qend, qstart)
                covered_pred_sections[pred].append(interval)
                covered_pred_sections_per_ref[(pred, ref)].append(interval)

    out.write("number of predicted plasmids: %i\n" % len(predicted_plasmids))
    for pred in predicted_plasmids:
        out.write("\t%s: %i nt\n" % (pred, len(predicted_sequences[pred])))
    out.write("\n")

    out.write("> predicted PLASMID covers <proportion> of reference plasmid\n")
    out.write("\t".join([""] + sorted(ref_plasmids)) + "\n")
    #out.write("\t".join([""] + sorted(ref_chromosomes)) + "\n")
    for pred in predicted_plasmids:
        out.write("\t".join([pred] + [str(num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])) for ref in sorted(ref_plasmids)]) + "\n")
        #out.write("\t".join([pred] + [str(num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])) for ref in sorted(ref_chromosomes)]) + "\n")
    out.write("\n")

    out.write("> reference plasmid covers <proportion> of predicted plasmid\n")
    out.write("\t".join([""] + sorted(predicted_sequences)) + "\n")
    for ref in sorted(ref_plasmids):
    #for ref in sorted(ref_chromosomes):
        out.write("\t".join([ref] + [str(num_covered_positions(covered_pred_sections_per_ref[(pred, ref)]) / len(predicted_sequences[pred])) for pred in sorted(predicted_sequences)]) + "\n")
    out.write("\n")

    out.write("> in total, how much of predicted plasmid is covered by reference plasmids\n")
    out.write("\t".join(["plasmid", "proportion"]) + "\n")
    for pred in sorted(predicted_sequences):
        out.write("%s\t%f\n" % (pred, num_covered_positions(covered_pred_sections[pred]) / len(predicted_sequences[pred])))
    out.write("\n")

    out.write("> in total, how much of reference plasmid is covered by predicted plasmids\n")
    out.write("\t".join(["plasmid", "proportion"]) + "\n")
    for ref in sorted(ref_plasmids):
    #for ref in sorted(ref_chromosomes):    
        out.write("%s\t%f" % (ref, num_covered_positions(covered_ref_sections[ref]) / len(reference_sequences[ref])) + "\n")
    out.write("\n")

    out.write("> pairs of predicted and reference plasmids with coverage >= %f in both directions\n" % threshold)
    empty = True
    for pred in predicted_plasmids:
        for ref in sorted(ref_plasmids):
        #for ref in sorted(ref_chromosomes):      
            proportion_ref = num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])
            proportion_potential = num_covered_positions(covered_pred_sections_per_ref[(pred, ref)]) / len(predicted_sequences[pred])
            if proportion_ref >= threshold and proportion_potential >= threshold:
                out.write("%s <-> %s (%f resp. %f)\n" % (ref, pred, proportion_ref, proportion_potential))
                empty = False
    if empty:
        out.write("none\n")
    out.write("\n")

    out.write("> summary scores (%s)\n" % name)
    # proportion of nucleotides of all reference plasmids covered by predicted plasmids
    len_to_cover = sum([len(reference_sequences[p]) for p in ref_plasmids])
    #len_to_cover = sum([len(reference_sequences[p]) for p in ref_chromosomes])
    len_covered = 0
    for ref in ref_plasmids:
    #for ref in ref_chromosomes:
        len_covered += num_covered_positions(covered_ref_sections[ref])
    recall = len_covered / len_to_cover if len_to_cover != 0 else 0
    out.write("recall: %f\n" % recall)

    # proportion of nucleotides of all predicted plasmids covered by reference plasmids
    len_to_cover = sum([len(predicted_sequences[pred]) for pred in predicted_sequences])
    len_covered = 0
    for pred in predicted_sequences:
        len_covered += num_covered_positions(covered_pred_sections[pred])
    precision = len_covered / len_to_cover if len_to_cover != 0 else 0
    out.write("precision: %f\n" % precision)

    # # weighted (by size of plasmid) average of best 1:1 correspondences between predicted and reference plasmids
    # total_len_refs = sum([len(reference_sequences[p]) for p in ref_plasmids])
    # correspondences = []
    # for ref in ref_plasmids:
    #     max_val = 0
    #     max_weight = 0
    #     for pred in predicted_plasmids:
    #         proportion_ref = num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])
    #         proportion_potential = num_covered_positions(covered_pred_sections_per_ref[(pred, ref)]) / len(predicted_sequences[pred])
    #         if (proportion_ref + proportion_potential) / 2 > max_val:
    #             max_val = (proportion_ref + proportion_potential) / 2
    #             max_weight = len(reference_sequences[ref]) / total_len_refs if total_len_refs != 0 else 0
    #     correspondences.append((max_val, max_weight))
    # score_correspondences = sum([value * weight for value, weight in correspondences])
    # out.write("score_correspondences: %f\n" % score_correspondences)

    f1_score = (2 * (recall * precision) / (recall + precision)) if recall + precision > 0 else 0
    out.write("f1_score: %f\n" % f1_score)
    out.write("\n")

    print_csv_row(name, sample_id, [recall, precision, f1_score], ";")


# analyse the matches between the contigs of predicted plasmids and the reference plasmids
def analyse_from_contigs(hits, predicted_plasmids, contig_sequences, reference_sequences, ref_chromosomes, ref_plasmids, threshold, name, sample_id, out):
    predicted_lengths = dict()
    for id in predicted_plasmids:
        predicted_lengths[id] = sum([len(contig_sequences[c]) for c in predicted_plasmids[id]])

    covered_ref_sections = dict([(ref, []) for ref in reference_sequences.keys()]) # of lists
    covered_pred_sections = dict([(pred, []) for pred in sorted(predicted_plasmids)])
    covered_ref_sections_per_pred = dict()
    covered_pred_sections_per_ref = dict()

    for id in sorted(predicted_plasmids):
        for ref in sorted(ref_plasmids):
            covered_pred_sections_per_ref[(id, ref)] = []
            covered_ref_sections_per_pred[(ref, id)] = []
            for contig in predicted_plasmids[id]:
                pred_ref_hits = hits.loc[hits.qseqid == contig].loc[hits.sseqid == ref]


                for index, row in pred_ref_hits.iterrows():
                    sstart = int(row[8])
                    send = int(row[9])
                    interval = (sstart, send) if sstart <= send else (send, sstart)
                    covered_ref_sections[ref].append(interval)
                    covered_ref_sections_per_pred[(ref, id)].append(interval)


                    qstart = int(row[6])
                    qend = int(row[7])
                    interval = (contig, qstart, qend) if qstart <= qend else (contig, qend, qstart)
                    covered_pred_sections[id].append(interval)
                    covered_pred_sections_per_ref[(id, ref)].append(interval)

    out.write("number of predicted plasmids: %i\n" % len(predicted_plasmids))
    for pred in predicted_plasmids:
        out.write("\t%s: %i nt\n" % (pred, predicted_lengths[pred]))
    out.write("\n")

    out.write("> predicted PLASMID covers <proportion> of reference plasmid\n")
    out.write("\t".join([""] + sorted(ref_plasmids)) + "\n")
    for id in sorted(predicted_plasmids):
        out.write("\t".join([str(id)] + [str(num_covered_positions(covered_ref_sections_per_pred[(ref, id)]) / len(reference_sequences[ref])) for ref in sorted(ref_plasmids)]) + "\n")
    out.write("\n")

    out.write("> reference plasmid covers <proportion> of predicted plasmid\n")
    out.write("\t".join([""] + list(map(str, sorted(predicted_plasmids)))) + "\n")
    for ref in sorted(ref_plasmids):
        covered_positions = dict()
        for id in sorted(predicted_plasmids):
            covered_positions[id] = sum([num_covered_positions([(qstart, qend) for contig, qstart, qend in covered_pred_sections_per_ref[(id, ref)] if contig == c]) for c in predicted_plasmids[id]])
        out.write("\t".join([ref] + [str(covered_positions[id] / predicted_lengths[id]) for id in sorted(predicted_plasmids)]) + "\n")
    out.write("\n")

    out.write("> in total, how much of predicted plasmid is covered by reference plasmids\n")
    out.write("\t".join(["plasmid", "proportion"]) + "\n")
    for id in sorted(predicted_plasmids):
        covered_positions = sum([num_covered_positions([(qstart, qend) for contig, qstart, qend in covered_pred_sections[id] if contig == c]) for c in predicted_plasmids[id]])
        out.write("%s\t%f\n" % (id, covered_positions / predicted_lengths[id]))
    out.write("\n")

    out.write("> in total, how much of reference plasmid is covered by predicted plasmids\n")
    out.write("\t".join(["plasmid", "proportion"]) + "\n")
    for ref in sorted(ref_plasmids):
        out.write("%s\t%f\n" % (ref, num_covered_positions(covered_ref_sections[ref]) / len(reference_sequences[ref])))
    out.write("\n")

    out.write("> pairs of predicted and reference plasmids with coverage >= %f in both directions\n" % threshold)
    empty = True
    for pred in predicted_plasmids:
        for ref in sorted(ref_plasmids):
            proportion_ref = num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])
            proportion_potential = sum([num_covered_positions([(qstart, qend) for contig, qstart, qend in covered_pred_sections_per_ref[(pred, ref)] if contig == c]) for c in predicted_plasmids[pred]]) / predicted_lengths[pred]
            if proportion_ref >= threshold and proportion_potential >= threshold:
                out.write("%s <-> %s (%f resp. %f)\n" % (ref, pred, proportion_ref, proportion_potential))
                empty = False
    if empty:
        out.write("none\n")
    out.write("\n")

    out.write("> summary scores (%s)\n" % name)
    # proportion of nucleotides of all reference plasmids covered by predicted plasmids
    len_to_cover = sum([len(reference_sequences[p]) for p in ref_plasmids])
    len_covered = 0
    for ref in ref_plasmids:
        len_covered += num_covered_positions(covered_ref_sections[ref])
    recall = len_covered / len_to_cover if len_to_cover != 0 else 0
    out.write("recall: %f\n" % recall)

    # proportion of nucleotides of all predicted plasmids covered by reference plasmids
    len_to_cover = sum([predicted_lengths[pred] for pred in predicted_plasmids])
    len_covered = 0
    for pred in predicted_plasmids:
        len_covered += sum([num_covered_positions([(qstart, qend) for contig, qstart, qend in covered_pred_sections[pred] if contig == c]) for c in predicted_plasmids[pred]])
    precision = len_covered / len_to_cover if len_to_cover != 0 else 0
    out.write("precision: %f\n" % precision)

    # # weighted (by size of plasmid) average of best 1:1 correspondences between predicted and reference plasmids
    # total_len_refs = sum([len(reference_sequences[p]) for p in ref_plasmids])
    # correspondences = []
    # for ref in ref_plasmids:
    #     max_val = 0
    #     max_weight = 0
    #     for pred in predicted_plasmids:
    #         proportion_ref = num_covered_positions(covered_ref_sections_per_pred[(ref, pred)]) / len(reference_sequences[ref])
    #         proportion_potential = sum([num_covered_positions([(qstart, qend) for contig, qstart, qend in covered_pred_sections_per_ref[(pred, ref)] if contig == c]) for c in predicted_plasmids[pred]]) / predicted_lengths[pred]
    #         if (proportion_ref + proportion_potential) / 2 > max_val:
    #             max_val = (proportion_ref + proportion_potential) / 2
    #             max_weight = len(reference_sequences[ref]) / total_len_refs if total_len_refs != 0 else 0
    #     correspondences.append((max_val, max_weight))
    # score_correspondences = sum([value * weight for value, weight in correspondences])
    # out.write("score_correspondences: %f\n" % score_correspondences)

    f1_score = (2 * (recall * precision) / (recall + precision)) if recall + precision > 0 else 0
    out.write("f1_score: %f\n" % f1_score)
    out.write("\n")

    print_csv_row(name, sample_id, [recall, precision, f1_score], ";")


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("sample_id", help = "")
    argparser.add_argument("chromosomes_fasta", help = "")
    argparser.add_argument("plasmids_fasta", help = "")
    argparser.add_argument("detailed_output", help = "")
    argparser.add_argument("-i", "--ilp_putative", default = "", help = "TODO") # name, queries, blast    
    argparser.add_argument("-p", "--greedy_putative", default = "", help = "TODO") # name, queries, blast
    argparser.add_argument("-q", "--greedy_questionable", default = "", help = "TODO")
    argparser.add_argument("-s", "--plasmidspades", default = "", help = "TODO")
    argparser.add_argument("-m", "--mob_recon", default = "", help = "TODO")
    argparser.add_argument("-t", "--threshold", type = float, default = 0.8, help = "TODO")
    argparser.add_argument("-c", "--greedy_contigs", action = "store_true", help = "TODO")
    argparser.add_argument("-d", "--ilp_contigs", action = "store_true", help = "TODO")
    args = argparser.parse_args()

    references = dict()
    chromosomes = set()
    plasmids = set()
    with open(args.chromosomes_fasta, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            references[record.id] = record.seq
            chromosomes.add(record.id)
    with open(args.plasmids_fasta, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            references[record.id] = record.seq
            plasmids.add(record.id)

    with open(args.detailed_output, "w") as out:

        out.write("##### general information #####\n")
        out.write("number of reference chromosomes: %i\n" % len(chromosomes))
        for chr in chromosomes:
            out.write("\t%s: %i nt\n" % (chr, len(references[chr])))
        out.write("number of reference plasmids: %i\n" % len(plasmids))
        for pla in plasmids:
            out.write("\t%s: %i nt\n" % (pla, len(references[pla])))
        out.write("\n")

        # ilp
        if args.ilp_putative != "":
            name, fasta_file, blast_file = args.ilp_putative.split(";")

            out.write("##### ILP (%s) #####\n" % name)
            ilp_putative = read_blast_output(blast_file)
            plasmid_seqs = dict()
            with open(fasta_file, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    plasmid_seqs[record.id] = record.seq

            if args.ilp_contigs:
                predicted_ids = set([name.split('_')[-1] for name in plasmid_seqs])
                predicted_plasmids = dict()
                for id in predicted_ids:
                    predicted_plasmids[id] = [name for name in plasmid_seqs if name.split('_')[-1] == id]
                analyse_from_contigs(ilp_putative, predicted_plasmids, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)
            else:
                analyse_greedy(ilp_putative, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)


        # greedy
        if args.greedy_putative != "":
            name, fasta_file, blast_file = args.greedy_putative.split(";")

            out.write("##### greedy putative (%s) #####\n" % name)
            greedy_putative = read_blast_output(blast_file)
            plasmid_seqs = dict()
            with open(fasta_file, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    plasmid_seqs[record.id] = record.seq

            if args.greedy_contigs:
                predicted_ids = set([name.split('_')[-1] for name in plasmid_seqs])
                predicted_plasmids = dict()
                for id in predicted_ids:
                    predicted_plasmids[id] = [name for name in plasmid_seqs if name.split('_')[-1] == id]
                analyse_from_contigs(greedy_putative, predicted_plasmids, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)
            else:
                analyse_greedy(greedy_putative, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)

        if args.greedy_questionable != "":
            name, fasta_file, blast_file = args.greedy_questionable.split(";")

            out.write("##### greedy questionable (%s) #####\n" % name)
            greedy_questionable = read_blast_output(blast_file)
            plasmid_seqs = dict()
            with open(fasta_file, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    plasmid_seqs[record.id] = record.seq

            if args.greedy_contigs:
                predicted_ids = set([name.split('_')[-1] for name in plasmid_seqs])
                predicted_plasmids = dict()
                for id in predicted_ids:
                    predicted_plasmids[id] = [name for name in plasmid_seqs if name.split('_')[-1] == id]
                analyse_from_contigs(greedy_questionable, predicted_plasmids, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)
            else:
                analyse_greedy(greedy_questionable, plasmid_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)


        # plasmidSPAdes
        if args.plasmidspades != "":
            name, fasta_file, blast_file = args.plasmidspades.split(";")

            out.write("##### plasmidSPAdes (%s) #####\n" % name)
            plasmidspades = read_blast_output(blast_file)
            contig_seqs = dict()
            with open(fasta_file, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    contig_seqs[record.id] = record.seq

            predicted_ids = set([name.split('_')[-1] for name in contig_seqs])
            predicted_plasmids = dict()
            for id in predicted_ids:
                predicted_plasmids[id] = [name for name in contig_seqs if name.split('_')[-1] == id]
            analyse_from_contigs(plasmidspades, predicted_plasmids, contig_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)


        # mob_recon
        if args.mob_recon != "":
            name, fasta_files, blast_files = args.mob_recon.split(";")
            fasta_files = fasta_files.split(',') if len(fasta_files) > 0 else []
            blast_files = blast_files.split(',') if len(blast_files) > 0 else []

            out.write("##### mob_recon (%s) #####\n" % name)
            mob_recon = pd.concat([read_blast_output(f) for f in blast_files], ignore_index = True)
            contig_seqs = dict()
            predicted_plasmids = dict()
            for file in fasta_files:
                plasmid_name = file[file.rfind("/") + 1 : file.rfind(".")]
                predicted_plasmids[plasmid_name] = []
                with open(file, "rU") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        contig_seqs[record.id] = record.seq
                        predicted_plasmids[plasmid_name].append(record.id)

            analyse_from_contigs(mob_recon, predicted_plasmids, contig_seqs, references, chromosomes, plasmids, args.threshold, name, args.sample_id, out)


main()