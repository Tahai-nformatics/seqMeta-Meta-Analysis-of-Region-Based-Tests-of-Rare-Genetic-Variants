#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
In collaboration with Chienyueh Lee

"""





from os import getcwd, path, makedirs
import sys
from optparse import OptionParser
import subprocess as sp
import gzip
import pandas as pd
from multiprocessing import Pool, cpu_count, Value
from glob import glob
from timeit import default_timer as timer
from datetime import datetime
import drmaa
import shlex
from pybgen import ParallelPyBGEN as PyBGEN
import errno

__version__ = '1.2.0'


def main():
    global config
    
    # Options
    usage = "\nStart a new analysis: %prog [options] <VCF file/-g GEN file/-b BGEN file> <Directory of SNP info files> <phenotype file> <chr number>\nResume previous work: %prog [options] <-r seqMeta_pipeline.log>\n\nVersion: "+__version__
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--cpu', dest='cpu', action='store', type='int', default=cpu_count(), help="How many jobs submitted simultaneously at one time [Default: %d]" % cpu_count())
    parser.add_option('-g', '--gen', dest='use_gen', action='store_true', default=False, help='Specify the input file is a GEN file instead of default VCF file [Default: False]')
    parser.add_option('-b', '--bgen', dest='use_bgen', action='store_true', default=False, help='Specify the input file is a BGEN file instead of default VCF file [Default: False]')
    parser.add_option('-p', '--prefix', dest='prefix', action='store', type='string', default='', help="Add a prefix tag in file names of all outputs [Default: '']")
    parser.add_option('-r', '--resume', dest='resume_log', action='store', type='string', default='', help="Give a seqMeta_pipeline.log file to resume previous work [Default: '']")
    parser.add_option('-u', '--unverify', dest='unverify', action='store_true', default=False, help='Skip to verify sample IDs and their order between VCF and phenotype files. (Refer to -v option) [Default: False]')
    parser.add_option('-v', '--verify', dest='verify_col_name', action='store', type='string', default='IID', help="Specify a column name of the phenotype file to verify sample IDs and their order between VCF and phenotype files. (e.g., '-v IID' indicates sample IDs in the IID column) [Default: 'IID']")
    parser.add_option('-w', '--workdir', dest='working_dir', action='store', type='string', default=getcwd(), help="The working directory for input/output files of seqMeta[Default: %s]" % getcwd())
    parser.add_option('-j', '--project', dest='project', action='store', type='string', default='', help="Project job tag [Default: '']")
    parser.add_option('--rscript_bin', dest='Rscript_bin_path', action='store', type='string', default='Rscript', help="Path of executable file for Rscript [Default: 'Rscript']")
    parser.add_option('--version', dest='version', action='store_true', default=False, help='Show the version of this script')

    (options, args) = parser.parse_args()
    
    if options.version:
        print(__version__)
        sys.exit(0)
    
    if not((len(options.resume_log) == 0 and len(args) == 4) or (len(options.resume_log) != 0 and len(args) == 0)):
        parser.print_help()
        sys.exit(1)
        
    if len(options.resume_log)==0: ## Start a new analysis
        ## Processing global config
        config.update({'chr': args[3],
                       'cpu': options.cpu,
                       'prefix': options.prefix,
                       'working_dir': options.working_dir,
                       'unverify': options.unverify,
                       'verify': options.verify_col_name,
                       'project': options.project,
                       'Rscript_bin_path': options.Rscript_bin_path,
                      })
        if len(config['prefix'])>0 : config['prefix'] += '.'
        if not path.isabs(config['working_dir']):
            config['working_dir'] = path.abspath(path.join('.', config['working_dir']))
        input_file = path.abspath(path.expanduser(args[0]))
        input_snp_info_dir = path.abspath(path.expanduser(args[1]))
        input_phenotype_file = path.abspath(path.expanduser(args[2]))
        
        if options.use_bgen:
            analysis_log_path = bgen2dosage(input_file, input_snp_info_dir, input_phenotype_file)
        else:
            if options.use_gen:
                gen_file = input_file
            else:
                ## Start to convert VCF to GEN
                gen_file = vcf2gen(input_file, input_phenotype_file)
            
            ## Start to process dosage files
            analysis_log_path = gen2dosage(gen_file, input_snp_info_dir, input_phenotype_file)
        
    else: ## Resume previous work
        ## Processing global config
        config.update({'cpu': options.cpu,
                       'Rscript_bin_path': options.Rscript_bin_path,
                      })
        analysis_log_path = options.resume_log
    
    ## Start to seqMeta analysis
    seqMeta(analysis_log_path)


def vcf2gen(input_vcf, input_phenotype_file):
    if not _is_path_exist(input_vcf, True): sys.exit(1)
    if not _is_path_exist(input_phenotype_file, True): sys.exit(1)
    
    _show_progress_message('Converting VCF to GEN files')
    
    ##　Create a directory for saving GEN files
    gen_file_dir = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_input', 'gen_file')
    _mkdir(gen_file_dir)
    gen_file_path = path.join(gen_file_dir, 'chr%s.gen.gz' % config['chr'])
    
    phenotype_df = pd.read_csv(input_phenotype_file, sep='\t')
    
    with gzip.open(gen_file_path, 'wt') as write_fp:
        input_vcf_basename = path.basename(input_vcf)    
        if path.splitext(input_vcf_basename)[1] == '.gz':
            read_fp = gzip.open(input_vcf, 'rt')
        else:
            read_fp = open(input_vcf, 'r')
        for line in read_fp:
            line = line.strip()
            if not line or line.startswith('##'): continue
            if line.startswith('#CHROM'):
                ## Verify the sample IDs and order
                if not config['unverify']:
                    vcf_header = line.split('\t')[9:]
                    phenotype_header = pd.read_csv(input_phenotype_file, delim_whitespace=True)[config['verify']].astype('str').values.tolist()
                    if vcf_header != phenotype_header:
                        sys.stderr.write("ERROR: Sample IDs between VCF and phenotype files are NOT identical. Please make sure IDs match before continuing with your analysis.\n")
                        sys.exit(1)
                
                continue
                
            ## line = ['22', '16050435', '22:16050435', 'T', 'C', '.', 'PASS', 'MAF=0.00043;R2=0.01842', 'GT:DS:GP', '0/0:0.001:0.999,0.001,0.000', ...]
            line = line.split('\t')
            line[1], line[2] = line[2], line[1]
            
            ## Check essential info
            if 'GP' not in line[8]:
                sys.stderr.write("Data not found: Estimated Posterior Probabilities for Genotypes (GP). Missing essential data in the VCF for converting a GEN format.\n")
                sys.exit(1)
            
            output_str1 = ' '.join(line[0:5])
            output_str2 = ' '.join(map(lambda x: ' '.join(x.split(':')[2].split(',')), line[9:]))
            ## Output gen format: 22 22:16050435 16050435 T C 0.001 0.999 0.000 0.010 0.980 0.010 ...
            write_fp.write('%s %s\n' % (output_str1, output_str2))
            
        read_fp.close()
        
    _show_progress_message('Completed VCF to GEN conversion')
    return gen_file_path
            

def gen2dosage(input_gen, input_snp_info_dir, input_phenotype_file):
    if not _is_path_exist(input_gen, True): sys.exit(1)
    
    _show_progress_message('Processing dosage files')
    
    dosage_file_pointers = {} ## {'GENE_ID': FILE_POINTER, ...}
    SNPinfo_pointers = {} ## {'GENE_ID': FILE_POINTER, ...}
    return_info = {} ## {'GENE': {'dosage_file_path': DOSAGE_FILE_PATH, 'SNPinfo_file_path': SNPINFO_FILE_PATH, 'snp_number': INT}}
    
    ##　Create a directory for saving dosage files
    dosage_file_dir = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_input', 'dosage_file')
    _mkdir(dosage_file_dir)
    
    ##　Create a directory for saving SNPinfo files
    SNPinfo_file_dir = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_input', 'SNPinfo')
    _mkdir(SNPinfo_file_dir)
            
    ## Load SNP info files
    snp_info_df = pd.read_csv(path.join(input_snp_info_dir, 'chr%s.snpinfo.gz' % config['chr']), compression='gzip', sep='\t')
    snp2gene = snp_info_df.groupby(['snpname']).aggregate(lambda x: tuple(x))['gene'].to_dict()
    gene_count = len(set(snp2gene.values()))

    ## The position table contains the maximum position of each gene getting from SNP info files. This table uses to determine those opened file pointers which can be closed at an appropriate time, to prevent exceeding the limitation of the maximum number of open files from ulimit 
    pos_df = snp_info_df['snpname'].str.split(':', 1, expand=True)[[1]].astype('int64').rename(columns={1: 'position'})
    snp_info_df = pd.concat([snp_info_df, pos_df], axis=1)
    max_pos_of_gene_df = snp_info_df.groupby(['gene']).position.agg('max').to_frame().reset_index()
    del(snp_info_df)
    
    input_gen_basename = path.basename(input_gen)    
    if path.splitext(input_gen_basename)[1] == '.gz':
        read_fp = gzip.open(input_gen, 'rt')
    else:
        read_fp = open(input_gen, 'r')

    rcd_count = 0
    for line in read_fp:
        ## 22 22:16050435 16050435 T C 1.000 0.000 0.000 0.996 0.004 0.000 ...
        line = line.strip().split(' ')
        snp_id = line[1]
        current_pos = int(line[2])
        ref = line[3]
        alt = line[4]
        
        ## The dosage information of genotypes starts at the 6th column with 3 columns as a sample
        genotype_dosage_list = line[5:]
        per_variant_list = ['%s%s>%s' % (snp_id, ref, alt)]
        for i in range(0, len(genotype_dosage_list), 3):
            each_sample = genotype_dosage_list[i:i+3]
            sample_dosage = float(each_sample[1]) + float(each_sample[2])*2
            per_variant_list.append(str(sample_dosage))
        
        try:
            for gene in snp2gene[snp_id]:
                if gene not in dosage_file_pointers:
                    rcd_count += 1
                    _show_progress_message("Creating dosage file {0}{1}.dosage.gz, gene {2} of {3} ({4:.2f}%) processed".format(config['prefix'], gene, rcd_count, gene_count, (100.0*rcd_count/gene_count)))
                    dosage_file_path = path.join(dosage_file_dir, '%s%s.dosage.gz' % (config['prefix'], gene))
                    SNPinfo_file_path = path.join(SNPinfo_file_dir, '%s%s.SNPinfo' % (config['prefix'], gene))
                    return_info[gene] = {'gene': gene, 'dosage_file_path': dosage_file_path, 'SNPinfo_file_path': SNPinfo_file_path, 'snp_number': 1}
                    
                    dosage_file_pointers[gene] = gzip.open(dosage_file_path, 'wb', compresslevel=9)
                    SNPinfo_pointers[gene] = open(SNPinfo_file_path, 'w')
                    SNPinfo_pointers[gene].write('snpname\tgene\n') ## Make header
                    
                    return_info[gene] = {'gene': gene, 'dosage_file_path': dosage_file_path, 'SNPinfo_file_path': SNPinfo_file_path, 'snp_number': 1}
                else:
                    return_info[gene]['snp_number'] += 1

                dosage_file_pointers[gene].write(('\t'.join(per_variant_list)+'\n').encode())
                SNPinfo_pointers[gene].write('%s%s>%s\t%s\n' % (snp_id, ref, alt, gene))
            
            ## Close file pointers according to maximum positions of genes which are less than the current position
            for need_close_gene in max_pos_of_gene_df.loc[max_pos_of_gene_df.position <= current_pos, 'gene']:
                if need_close_gene in dosage_file_pointers:
                    dosage_file_pointers[need_close_gene].close()
                    SNPinfo_pointers[need_close_gene].close()
                    del(dosage_file_pointers[need_close_gene])
                    del(SNPinfo_pointers[need_close_gene])
            
        except KeyError:
            continue ## Ignore SNPs w/o in the snp_list
    
    if rcd_count < gene_count:
        _show_progress_message("Skipped {0} genes without SNPs, gene {1} of {2} (100.00%) processed".format((gene_count-rcd_count), gene_count, gene_count))

    ## Close remaining file pointers
    for fp in dosage_file_pointers.values():
        fp.close()
    for fp in SNPinfo_pointers.values():
        fp.close()

    read_fp.close()
    
    ## Save analysis log
    select_col_by_order = ['chr', 'gene', 'snp_number', 'dosage_file_path', 'SNPinfo_file_path', 'phenotype_file_path', 'prefix', 'working_dir', 'project', 'job_id', 'process_status']
    path_df = pd.DataFrame([return_info[gene] for gene in return_info.keys()]).sort_values(by=['snp_number'])
    path_df['phenotype_file_path'] = input_phenotype_file
    path_df['chr'] = config['chr']
    path_df['prefix'] = config['prefix']
    path_df['working_dir'] = config['working_dir']
    path_df['project'] = config['project']
    path_df['job_id'] = None
    path_df['process_status']= 'waiting'
    analysis_log = path.join(config['working_dir'], 'chr'+config['chr'], '%schr%s_before_seqMeta_analysis.log' % (config['prefix'], config['chr']))
    path_df[select_col_by_order].to_csv(analysis_log, sep='\t', index=False)
    
    return analysis_log


def _get_variants_by_region(input_bgen, gene_region_list, gene_count):
    global processed_gene_num
    ## gene_region_list = ['GENE_ID', 'chr', pos_start, pos_end]
    gene_id = gene_region_list[0]
    gene_chr = str(gene_region_list[1]).replace('chr', '')
    gene_pos_start = gene_region_list[2]
    gene_pos_end = gene_region_list[3]
    snp_number = 0
    dosage_list = []   ## 
    SNPinfo_list = ["snpname\tgene"]
    bgen = PyBGEN(input_bgen)
    
    for v,d in bgen.iter_variants_in_region(gene_chr, gene_pos_start, gene_pos_end):
        snp_number += 1
        snpname = '%s:%d%s>%s' % (v.chrom, v.pos, v.a1, v.a2)
        dosage_element_list = [snpname]
        
        dosage_element_list.extend(list(map(str, d)))
        dosage_list.append("\t".join(dosage_element_list))
        SNPinfo_list.append("%s\t%s" % (snpname, gene_id))
    
    if snp_number > 0:
        ## Count how many genes are processed
        with processed_gene_num.get_lock():
            processed_gene_num.value += 1
        ##　Create a directory for saving dosage files
        dosage_file_dir = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_input', 'dosage_file')
        _mkdir(dosage_file_dir)
        ##　Create a directory for saving SNPinfo files
        SNPinfo_file_dir = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_input', 'SNPinfo')
        _mkdir(SNPinfo_file_dir)
        ## File paths
        dosage_file_path = path.join(dosage_file_dir, '%s%s.dosage.gz' % (config['prefix'], gene_id))
        SNPinfo_file_path = path.join(SNPinfo_file_dir, '%s%s.SNPinfo' % (config['prefix'], gene_id))
        
        _show_progress_message("Creating dosage file {0}{1}.dosage.gz, gene {2} of {3} ({4:.2f}%) processed".format(config['prefix'], gene_id, processed_gene_num.value, gene_count, (100.0*processed_gene_num.value/gene_count)))
        with gzip.open(dosage_file_path, 'wb', compresslevel=9) as dosage_fp:
            dosage_fp.write("\n".join(dosage_list).encode())
        with open(SNPinfo_file_path, 'w') as SNPinfo_fp:
            SNPinfo_fp.write("\n".join(SNPinfo_list))
        
        return {'gene': gene_id, 'dosage_file_path': dosage_file_path, 'SNPinfo_file_path': SNPinfo_file_path, 'snp_number': snp_number}
        
    else: return {}


def bgen2dosage(input_bgen, input_snp_info_dir, input_phenotype_file):
    if not _is_path_exist(input_bgen, True): sys.exit(1)
    
    _show_progress_message('Processing dosage files')
    
    global processed_gene_num
    ## Load SNP info files and process
    snp_info_df = pd.read_csv(path.join(input_snp_info_dir, 'chr%s.snpinfo.gz' % config['chr']), compression='gzip', sep='\t')
    pos_df = snp_info_df['snpname'].str.split(':', 1, expand=True)[[0,1]].rename(columns={0: 'chr', 1: 'position'}).astype({'chr':'str', 'position':'int64'})
    snp_info_df = pd.concat([snp_info_df, pos_df], axis=1, sort=False)
    gene_region_df = snp_info_df.groupby(['gene']).agg({'chr':'first', 'position':['min', 'max']}).reset_index()
    del(snp_info_df)
    
    ## Prepare data for multi threads processing
    ## task_list = [['GENE_ID', 'chr', pos_start, pos_end], ...]
    task_list = gene_region_df[['gene', 'chr', 'position']].values.tolist()
    gene_count = len(task_list)
    processed_gene_num = Value('i', 0)
    pool = Pool(config['cpu'])
    res_pool_list = [pool.apply_async(_get_variants_by_region, (input_bgen,x,gene_count,)) for x in task_list]
    ## Collect return
    updated_dosage_file_list = [r.get() for r in res_pool_list]
    
    if processed_gene_num.value < gene_count:
        _show_progress_message("Skipped {0} genes without SNPs, gene {1} of {2} (100.00%) processed".format((gene_count-processed_gene_num.value), gene_count, gene_count))
    del(processed_gene_num)
    
    ## Save analysis log
    select_col_by_order = ['chr', 'gene', 'snp_number', 'dosage_file_path', 'SNPinfo_file_path', 'phenotype_file_path', 'prefix', 'working_dir', 'project', 'job_id', 'process_status']
    path_df = pd.DataFrame([y for y in updated_dosage_file_list if bool(y)]).sort_values(by=['snp_number'])
    path_df['phenotype_file_path'] = input_phenotype_file
    path_df['chr'] = config['chr']
    path_df['prefix'] = config['prefix']
    path_df['working_dir'] = config['working_dir']
    path_df['project'] = config['project']
    path_df['job_id'] = None
    path_df['process_status']= 'waiting'
    analysis_log = path.join(config['working_dir'], 'chr'+config['chr'], '%schr%s_before_seqMeta_analysis.log' % (config['prefix'], config['chr']))
    path_df[select_col_by_order].to_csv(analysis_log, sep='\t', index=False)
    
    return analysis_log


def _run_seqMeta(input_list, total_jobs):
    ## The schema of input_list:
    ## 0:'chr', 1:'gene', 2:'snp_number', 3:'dosage_file_path', 4:'SNPinfo_file_path', 5:'phenotype_file_path', 6:'prefix', 7:'working_dir', 8:'project', 9:'job_id', 10:'process_status'
    return_dict = {'chr':input_list[0], 'gene':input_list[1], 'snp_number':input_list[2], 'dosage_file_path':input_list[3], 'SNPinfo_file_path':input_list[4], 'phenotype_file_path':input_list[5], 'prefix':input_list[6], 'working_dir':input_list[7], 'project':input_list[8], 'job_id':input_list[9], 'process_status':'running'}
    seqMeta_output = path.join(return_dict['working_dir'], 'chr'+str(return_dict['chr']), 'seqMeta_output')
    job_log = path.join(return_dict['working_dir'], 'chr'+str(return_dict['chr']), 'logs')
    global processed_job_num
    
    returned_job_dict = submitJob(gene=return_dict['gene'],
                                  dosage=return_dict['dosage_file_path'],
                                  snpinfo=return_dict['SNPinfo_file_path'],
                                  pheno=return_dict['phenotype_file_path'],
                                  count=return_dict['snp_number'],
                                  prefix=return_dict['prefix'],
                                  seqMeta_output=seqMeta_output,
                                  job_log=job_log,
                                  project=return_dict['project'],
                                  chr=return_dict['chr'],)
    
    ## Count how many jobs are processed
    with processed_job_num.get_lock():
        processed_job_num.value += 1
    return_dict['job_id'] = returned_job_dict['jobId']
    submission_output_message = "JOB-ID: {} ".format(returned_job_dict['jobId'])
    if returned_job_dict['hasExited']:
        if returned_job_dict['exitStatus'] == 0:
            return_dict['process_status'] = 'Successfully completed'
            submission_output_message += 'is successfully completed'
        else:
            submission_output_message += 'is failed due to '
            if returned_job_dict['wasAborted']:
                return_dict['process_status'] = 'Job is aborted'
                submission_output_message += 'job aborted'
            elif returned_job_dict['hasCoreDump']:
                return_dict['process_status'] = 'Job has core dump'
                submission_output_message += 'core dump'
            elif returned_job_dict['hasSignal']:
                return_dict['process_status'] = "Job is terminated by signal {}".format(returned_job_dict['terminatedSignal'])
                submission_output_message += "the job is terminated by signal {}".format(returned_job_dict['terminatedSignal'])
            else:
                return_dict['process_status'] = 'Cannot allocate memory'
                submission_output_message += 'the job cannot allocate memory'
    else:
        return_dict['process_status'] = 'Job is not exited'
        submission_output_message += 'the job is not exited'
    
    submission_output_message += ", job {0} of {1} ({2:.2f}%) processed".format(processed_job_num.value, total_jobs, (100.0*processed_job_num.value/total_jobs))
    _show_progress_message(submission_output_message)
    
    return return_dict


def seqMeta(analysis_log):
    if not _is_path_exist(analysis_log, True): sys.exit(1)
    
    _show_progress_message('Performing seqMeta analysis')
    
    dosage_file_df = pd.read_csv(analysis_log, sep='\t', engine='python').fillna('')
    ## Get values from analysis_log instead of config because of the resume function
    chr = str(dosage_file_df['chr'][0])
    prefix = dosage_file_df['prefix'][0]
    global processed_job_num
    
    ## Select genes at least comprising 3 SNPs as seqMeta inputs. Order the columns and convert to a list for input
    select_col_by_order = ['chr', 'gene', 'snp_number', 'dosage_file_path', 'SNPinfo_file_path', 'phenotype_file_path', 'prefix', 'working_dir', 'project', 'job_id', 'process_status']
    task_list = dosage_file_df[(dosage_file_df['snp_number']>2) & (dosage_file_df['process_status']!='Successfully completed')][select_col_by_order].values.tolist()
    total_jobs = len(task_list)
    
    ## Run seqMeta in parallel
    processed_job_num = Value('i', 0)
    pool = Pool(config['cpu'])
    res_pool_list = [pool.apply_async(_run_seqMeta, (x,total_jobs,)) for x in task_list]
    ## Collect return
    updated_dosage_file_list = [r.get() for r in res_pool_list]
    
    ## Create after_seqMeta_analysis.log
    after_seqMeta_analysis_log = path.join(path.dirname(analysis_log), '%schr%s_after_seqMeta_analysis_%s.log') % (prefix, chr, datetime.now().strftime("%m%d%Y%H%M%S")) ## Save as a new file
    
    updated_dosage_file_df = dosage_file_df[dosage_file_df['snp_number']<=2][:]
    updated_dosage_file_df['process_status'].replace('waiting', 'skipped', inplace=True)
    
    if len(updated_dosage_file_list)>0:
        tmp_df = pd.DataFrame(updated_dosage_file_list)
        updated_dosage_file_df = pd.concat([updated_dosage_file_df, tmp_df], sort=False).sort_values(by=['snp_number'])
        updated_dosage_file_df.to_csv(after_seqMeta_analysis_log, sep='\t', index=False)


def _is_path_exist(dir, error_msg=False):
    if path.exists(dir):
        return True
    else:
        if error_msg: sys.stderr.write(dir+': No such file or directory\n')
        return False


def _mkdir(dir):
    if not path.exists(dir):
        try:
            makedirs(dir)
        except OSError as err:
            if err.errno == errno.EEXIST: pass
            else: raise


def _show_progress_message(msg):
    sys.stdout.write("[%s] %s\n" % (datetime.now().strftime("%m-%d-%Y %H:%M:%S"), msg))
    sys.stdout.flush()


def submitJob(gene, dosage, snpinfo, pheno, count, prefix, seqMeta_output, job_log, project, chr):
    """Submit a job on the queue system using drmaa

    :param gene: Gene symbol
    :param dosage: path to dosage file
    :param snpinfo: path to snpinfo file
    :param pheno: path to phenotypes file
    :param count: count of markers, for memory estimation
    :param prefix: Prefix tag for seqMeta output files
    :param seqMeta_output: seqMeta outputs dir
    :param job_log: Job log dir
    :param project: Project job tag
    :returns: dict
    """
    return_dict = {}
    # Memory estimation
    mem = 2048 if count < 2048 else int(round(count*count*8*10*1.5/1024/1024, 0))

    with drmaa.Session() as s:
        # seqMeta_output = path.join(config['working_dir'], 'chr'+config['chr'], 'seqMeta_output/')
        # job_log = path.join(config['working_dir'], 'chr'+config['chr'], 'logs')
        _mkdir(seqMeta_output)
        _mkdir(job_log)
        prj = '-' + project if bool(project) else ''

        cmd_args = '{RPATH} {GENE} {D} {S} {P} "{PREFIX}" {OUT}'.format(
                     RPATH  = config['rscript_path'],
                     GENE   = gene,
                     D      = dosage,
                     S      = snpinfo,
                     P      = pheno,
                     PREFIX = prefix,
                     OUT    = seqMeta_output
                     )

        jt = s.createJobTemplate()
        jt.outputPath = ':' + job_log + '/' # log path (dir)
        jt.joinFiles=True
        jt.jobName = 'SeqMeta{}-chr{}-{}'.format(prj, chr, gene)
        jt.remoteCommand = config['Rscript_bin_path']
        jt.args = shlex.split(cmd_args)

        # Additional arguments to submit command
        if 'SGE' in s.drmsInfo or 'OGS' in s.drmsInfo:
            jt.nativeSpecification = '-l h_vmem=%dM' % (mem)
        elif 'LSF' in s.drmsInfo:
            jt.nativeSpecification = '-M %d -R "span[hosts=1]"' % (mem)

        jobid = s.runJob(jt)
        _show_progress_message("DRMAA using {}. Your job has been submitted with JOB-ID: {}".format(s.drmsInfo, jobid))

        retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        return_dict = {'jobId': jobid, 'exitStatus': retval.exitStatus, 'wasAborted': retval.wasAborted, 'hasCoreDump': retval.hasCoreDump, 'hasExited': retval.hasExited, 'hasSignal': retval.hasSignal, 'terminatedSignal': retval.terminatedSignal}
    
    return return_dict


if __name__ == '__main__':
    config = {}
    config['rscript_path'] = path.join(path.dirname(path.realpath(__file__)), 'run_seqmeta.R')
    if not _is_path_exist(config['rscript_path'], True): sys.exit(1)

    start_time = timer()
    main()
    print ('[%s] All done (check log files for errors if any). %s sec. elapsed' % (datetime.now().strftime("%m-%d-%Y %H:%M:%S"), str(timer() - start_time)))
