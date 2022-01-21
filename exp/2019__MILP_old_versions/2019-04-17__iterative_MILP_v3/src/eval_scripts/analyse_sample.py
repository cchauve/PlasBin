#!/usr/bin/python

# MILP strategy
#

import os
import subprocess


proj_dir = '/home/amane/projects/rrg-chauvec/wg-anoph/Plasmids-Optimization'
compare_script = proj_dir + '/exp/2019-03-22__modified_iterative_MILP/src/eval_scripts/compare_plasmids.py'
extract_script = proj_dir + '/data/2018-05-17__MOB-suite_benchmark/extract.sh'
sample_script = proj_dir + '/data/2018-05-23__MOB-suite_benchmark_reads/sample.sh'

#load_modules = 'module load gcc/5.4.0 blast+/2.6.0 boost/1.60.0 sickle/1.33 samtools/1.5 bowtie2/2.3.3.1 racon/20170719 perl/5.22.4 java/1.8.0_121 python/3.5.4; source $HOME/py3.5.4/bin/activate;'
#unload_modules = 'source deactivate; module unload python/3.5.4 java/1.8.0_121 perl/5.22.4 boost/1.60.0 blast+/2.6.0 gcc/5.4.0;'
    
def evaluate(sample_id, out_dir, results):
    p = subprocess.Popen('bash %s chr %s' % (sample_script, sample_id), stdout = subprocess.PIPE, shell = True)
    output, _ = p.communicate()
    p.wait()
    chromosomes = output.rstrip().decode().split(',')    

    p = subprocess.Popen('bash %s pla %s' % (sample_script, sample_id), stdout = subprocess.PIPE, shell = True)
    output, _ = p.communicate()
    p.wait()
    plasmids = output.rstrip().decode().split(',')    
    
    eval_dir = out_dir + '/eval'
    subprocess.call('mkdir %s' % eval_dir, shell = True)
    chr_references_fasta = eval_dir + '/chromosomes.fasta'
    pla_references_fasta = eval_dir + '/plasmids.fasta'
    subprocess.call('rm -f %s %s' % (chr_references_fasta, pla_references_fasta), shell = True)
    for acc in chromosomes:
        file = '%s/%s.fasta' % (out_dir, acc)
        subprocess.call('bash %s %s %s; cat %s >> %s; rm %s' % (extract_script, acc, file, file, chr_references_fasta, file), shell = True)
    for acc in plasmids:
        file = '%s/%s.fasta' % (out_dir, acc)
        subprocess.call('bash %s %s %s; cat %s >> %s; rm %s' % (extract_script, acc, file, file, pla_references_fasta, file), shell = True)

    blast_db = eval_dir + '/plasmid_references'
    subprocess.call('makeblastdb -in %s -dbtype nucl -out %s; ' % (pla_references_fasta, blast_db), shell = True)
    
    quast_labels = []
    quast_files = []
    
    
    ## evaluate MILP-based plasmids
    # first, map plasmid genes to Unicycler assembly and filter the mapping
    print('Performing evaluations with MILP strategy...')
    subprocess.call('mkdir %s/MILP' % eval_dir, shell = True)
    
    # MILP / partial training                     
    MILP_putative_queries = out_dir + '/MILP/putative_plasmids.fasta'
    MILP_putative_mapping = eval_dir + '/MILP/MILP_putative_map.csv'
    #greedy_questionable_queries = out_dir + '/greedy/plasmids/greedy/questionable_plasmids.fasta'
    #greedy_questionable_mapping = eval_dir + '/greedy/greedy_questionable_map.csv'
    subprocess.call( 'blastn -task megablast -query %s -db %s -out %s -outfmt 6; ' % (MILP_putative_queries, blast_db, MILP_putative_mapping) \
                    #+ 'blastn -task megablast -query %s -db %s -out %s -outfmt 6; ' % (greedy_questionable_queries, blast_db, greedy_questionable_mapping) \
                    + 'python %s %s %s %s %s/MILP/MILP_eval.csv -i "MILP_putative;%s;%s" >> %s; ' % (compare_script, sample_id, chr_references_fasta, pla_references_fasta, eval_dir, MILP_putative_queries, MILP_putative_mapping, results) , shell = True)

if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('sample_id', help = 'sample id of MOB-suite sample')
    argparser.add_argument('out_dir', help = 'output directory')
    argparser.add_argument('results', help = 'CSV file of scores')
    args = argparser.parse_args()
    
    evaluate(args.sample_id, args.out_dir, args.results)

    print('Finished analysis of sample %s.' % args.sample_id)
