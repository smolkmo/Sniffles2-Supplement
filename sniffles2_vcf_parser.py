#!/usr/bin/env python3

# Script parse sniffles2 vcf files and extract useful information
# USAGE:  cat sv_file.vcf | python3 sniffles2_vcf_parser.py


# important imports
import sys
import argparse
import random


# Version
def the_version():
    print(f'sniffles2_vcf_parser.py\nversion 0.1')

# log output
def logger(text, logtype=""):
    if logtype == "warn":
        sys.stderr.write("[WARNING] %s\n" % (text))
    elif logtype == "dev":
        sys.stderr.write("[DEV] %s\n" % (text))
    elif logtype == "error":
        sys.stderr.write("[ERROR] %s\n" % (text))
        sys.exit("Program stopped.")
    else:
        sys.stdout.write("[INFO] %s\n" % (text))

# vcf header class
class vcf_header(object):
    # init
    def __init__(self, input_line):
        # vcf files have 9 fields + the indiviuals data
        # split and analysis
        self.VCF_MANDATORY_FIELDS = 9 # not counting individuals genotype
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) < self.VCF_MANDATORY_FIELDS   # 9 VCF fields + genotypes
        if not self.ERROR:
            #extract indiviual(s): 1 -> n
            self.SAMPLES = tab_sep_fields[self.VCF_MANDATORY_FIELDS:]
        else:
            self.SAMPLES = ""

    def get_sample_name(self, svector, gt_vector=[], id2fam={}):
        all_names = []
        # extract names
        if len(gt_vector) > 0:
            for idx in range(0, len(self.SAMPLES)):
                supp=svector[idx]
                sname=self.SAMPLES[idx]
                gt_smp=gt_vector[idx]
                if supp == "1":
                    fmemb = ""
                    if sname in id2fam:
                        fmemb = id2fam[sname]
                    all_names.append(f'{sname}({gt_smp}){fmemb}')
        else:
            for sname, supp in zip(self.SAMPLES, list(svector)):
                if supp == "1":
                    all_names.append(sname)
        # make names vector into string
        if len(all_names) == 1:
                return all_names[0]
        else:
            return ", ".join(all_names)

    def get_sample_name_famtrio(self, gt_vector, id2fam):
        all_names = []
        fam_occurrence = {}
        # extract names
        if len(gt_vector) > 0:
            for sname, gt_smp in zip(self.SAMPLES, gt_vector):
                fmemb = ""
                if sname in id2fam:
                    fmemb = id2fam[sname]
                    fam_occurrence[fmemb] = gt_smp
                all_names.append(f'{sname}')
            # make names vector into string
            return [", ".join(all_names), fam_occurrence]
        else:
            # no info
            return ["", {}]


# vcf line class
class vcf_line(object):
    # init
    def __init__(self, input_line):
        # save the original
        self.svline = input_line
        # vcf files have 9 fields + the genotype data from the analyzed genome
        # we expect a single sample in this VCF file, for population there is another class
        # split and analysis
        self.VCF_MANDATORY_FIELDS = 9
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) < self.VCF_MANDATORY_FIELDS   # 9 VCF fields + genotype
        if not self.ERROR:
            #extract mandatory fileds
            [self.CHROM, POS, self.ID, REF, ALT, QUAL, self.FILTER, self.INFO, FORMAT] = tab_sep_fields[:self.VCF_MANDATORY_FIELDS]
            self.POS = int(POS)
            self.AF = -1
            self.GENOTYPE = ""
            self.get_genotype(tab_sep_fields, FORMAT)
            self.get_parsed_info()
            self.AF = self.DV/(self.DR+self.DV) if (self.AF == -1 and self.DR+self.DV > 0) else self.AF
            self.TRA = "" if self.SVTYPE != "BND" else ALT
            self.SVLEN = int(self.SVLEN) if self.SVTYPE != "BND" else self.SVLEN

    def get_genotype(self, tab_sep_fields, FORMAT):
        # FORMAT field and genotype extraction
        ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">
        ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">
        ## AF : Allele frequency. AF can be in FORMAT or INFO :(
        ## DP : Depth
        extract_genotype_info = ["GT", "DR", "DV", "AF", "DP"]
        self.DR = 0
        self.DV = 0
        self.DP = 0
        if ":" in FORMAT:
            # indiviuals in the vcf file
            gt_format_list = FORMAT.split(":")
            gt_list = tab_sep_fields[self.VCF_MANDATORY_FIELDS].split(":")
            for frmt_field, sv_gt in zip(gt_format_list, gt_list):
                self.GENOTYPE = sv_gt if frmt_field == "GT" else self.GENOTYPE
                self.AF = float(sv_gt) if frmt_field == "AF" else self.AF
                self.DR = int(sv_gt) if frmt_field == "DR" else self.DR
                self.DV = int(sv_gt) if frmt_field == "DV" else self.DV
                self.DP = int(sv_gt) if frmt_field == "DP" else self.DP
        else:
            self.GENOTYPE = tab_sep_fields[self.VCF_MANDATORY_FIELDS]

    def get_parsed_info(self):
        # INFO field extraction
        # interested in SVTYPE=<STRING>;SVLEN=<SIGNED INT>;...;RE=<INT>
        ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
        ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
        ##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
        ##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">
        # SUPP=2;SUPP_VEC=11;SVLEN=293;SVTYPE=INS;SVMETHOD=SURVIVOR1.0.7;CHR2=chr1;END=136556;CIPOS=0,380;CIEND=0,527;STRANDS=+-
        extract_info = ["SVLEN", "SVTYPE", "SUPP_VEC", "AF", "AC", "SUPPORT", "RE"]
        self.REF_COUNT = 0
        self.ALT_COUNT = 0
        self.SUPPORT = 0
        self.SVLEN = "."
        self.SVTYPE = ""
        self.SUPP_VEC = ""
        if ";" in self.INFO:
            for each_info in self.INFO.split(";"):
                if "=" in each_info:
                    try:
                        [key, val] = each_info.split("=")
                    except:
                        logger(each_info, "warn")
                    self.SUPPORT = int(val) if (key == "SUPPORT" or key == "RE") else self.SUPPORT
                    self.SVLEN = val if key == "SVLEN" else self.SVLEN
                    self.SVTYPE = val if key == "SVTYPE" else self.SVTYPE
                    self.SUPP_VEC = val if key == "SUPP_VEC" else self.SUPP_VEC
                    self.AF = float(val) if key == "AF" else self.AF
                    if key == "AC" and "," in val:
                        [tmpREF, tmpALT] = val.split(",")
                        self.REF_COUNT = int(tmpREF)
                        self.ALT_COUNT = int(tmpALT)
                        self.DP = self.REF_COUNT + self.ALT_COUNT
                        self.AF = self.ALT_COUNT / (self.REF_COUNT + self.ALT_COUNT)

    def sv_print(self):
        print(self.svline, end="\n")


# vcf line class for population file
class vcf_line_pop(object):
    # init
    def __init__(self, input_line):
        # save the original
        self.svline = input_line
        # vcf files have 9 fields + the indiviuals data
        # split and analysis
        self.VCF_MANDATORY_FIELDS = 9
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) < self.VCF_MANDATORY_FIELDS   # 9 VCF fields + genotypes
        if not self.ERROR:
            #extract mandatory fileds
            [self.CHROM, POS, self.ID, REF, ALT, QUAL, FILTER, self.INFO, FORMAT] = tab_sep_fields[:self.VCF_MANDATORY_FIELDS]
            self.POS = int(POS)
            self.GENOTYPE = ""
            self.get_genotype(tab_sep_fields, FORMAT)
            self.get_parsed_info()
            self.TRA = "" if self.SVTYPE != "BND" else ALT
            self.SVLEN = int(self.SVLEN) if self.SVTYPE != "BND" else 1
            self.END = int(self.END) if self.SVTYPE != "BND" else self.POS+1

    def get_genotype(self, tab_sep_fields, FORMAT):
        # FORMAT field and genotype extraction
        # Always used FORMAT fields
        # use FORMAT duhh!!!
        self.all_GT = []
        if ":" in FORMAT:
            # indiviuals in the vcf file
            split_format = FORMAT.split(":")
            for each_gt in tab_sep_fields[self.VCF_MANDATORY_FIELDS:]:
                split_gt = each_gt.split(":")
                for gtf,sv_gt in zip(split_format, split_gt):
                    if gtf == "GT":
                        # cant avoid ./. :(
                        self.all_GT.append(sv_gt)
                    else:
                        pass
            use_gt = ""
            for gt in self.all_GT:
                if use_gt == "":
                    # first case
                    use_gt = gt
                elif use_gt == "./." and gt != "./.":
                    # avoid not determined
                    use_gt = gt
                elif use_gt == "0/0" and (gt != "0/0" and gt != "./."):
                    # not ref homozygous
                    use_gt = gt
                elif use_gt == "0/1" and (gt == "1/1" and gt != "0/0" and gt != "./."):
                    # not ref homozygous
                    use_gt = gt
                else:
                    # 0/0 -> 0/0 and 0/1 -> 0/1
                    pass
            self.GENOTYPE = use_gt

        else:
            self.GENOTYPE = tab_sep_fields[self.VCF_MANDATORY_FIELDS]

    def get_parsed_info(self):
        # INFO field extraction
        # interested in SVTYPE=<STRING>;SVLEN=<SIGNED INT>;...;RE=<INT>
        ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
        ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
        ##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
        ##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">
        # SUPP=2;SUPP_VEC=11;SVLEN=293;SVTYPE=INS;SVMETHOD=SURVIVOR1.0.7;CHR2=chr1;END=136556;CIPOS=0,380;CIEND=0,527;STRANDS=+-
        extract_info = ["SVLEN", "SVTYPE", "END", "SUPP", "SUPP_VEC"]
        self.SVLEN = "."
        self.SVTYPE = ""
        self.SUPP = ""
        self.SUPP_VEC = ""
        self.END = ""
        if ";" in self.INFO:
            for each_info in self.INFO.split(";"):
                if "=" in each_info:
                    try:
                        [key, val] = each_info.split("=")
                    except:
                        logger(each_info, "warn")
                    self.SVLEN = val if key == "SVLEN" else self.SVLEN
                    self.SVTYPE = val if key == "SVTYPE" else self.SVTYPE
                    self.SUPP = val if key == "SUPP" else self.SUPP
                    self.SUPP_VEC = val if key == "SUPP_VEC" else self.SUPP_VEC
                    self.END = val if key == "END" else self.END

    def sv_print(self):
        print(self.svline, end="\n")


def get_genotype_sv(args):
    if args.only_genotype:
        print(f'#GENOTYPE\tCOUNT', end="\n")
    else:
        print(f'#CHROM\tPOSITION\tGENOTYPE\tSVTYPE\tSVLEN\tCOVERAGE\tAF\tREF\tALT', end="\n")
    gt_dict = {}
    gt_filter = [] if args.filer_gt == "" else args.filer_gt.split(",")
    for line in sys.stdin:
        # skip comments
        if not line.startswith("#"):
            vcf_entry = vcf_line(line.rstrip("\n"))
            gt = vcf_entry.GENOTYPE
            if gt not in gt_filter:
                if args.only_genotype:
                    if gt not in gt_dict:
                        gt_dict[gt] = 1
                    else:
                        gt_dict[gt] += 1
                else:
                    if (vcf_entry.DR+vcf_entry.DV) >= args.minsupp_read:
                        af = "{:0.3f}".format(vcf_entry.AF) if vcf_entry.AF != -1 else "NA"
                        if vcf_entry.SVTYPE == "BND":
                            print(f'{vcf_entry.CHROM}\t{vcf_entry.POS}\t{vcf_entry.GENOTYPE}\t{vcf_entry.SVTYPE}\t{vcf_entry.TRA}\t{vcf_entry.DR+vcf_entry.DV}\t{af}\t{vcf_entry.DR}\t{vcf_entry.DV}', end="\n")
                        else:
                            print(f'{vcf_entry.CHROM}\t{vcf_entry.POS}\t{vcf_entry.GENOTYPE}\t{vcf_entry.SVTYPE}\t{vcf_entry.SVLEN}\t{vcf_entry.DR+vcf_entry.DV}\t{af}\t{vcf_entry.DR}\t{vcf_entry.DV}', end="\n")
    if args.only_genotype:
        gt_total = 0
        for gt in gt_dict:
            print(f'{gt}\t{gt_dict[gt]}', end="\n")
            gt_total += gt_dict[gt]
        print(f'TOTAL\t{gt_total}', end="\n")


def population_parse(args):

    def print_out(vcf_entry, nSupport, support_ids_notRef):
        if vcf_entry.SVTYPE == "BND":
            print(f'{vcf_entry.CHROM}\t{vcf_entry.POS}\t{vcf_entry.END}\t{vcf_entry.GENOTYPE}\t{vcf_entry.SVTYPE}\t{vcf_entry.SVLEN}\t{nSupport}\t{support_ids_notRef}', end="\n")
        else:
            if abs(vcf_entry.SVLEN) > args.minsize:
                print(f'{vcf_entry.CHROM}\t{vcf_entry.POS}\t{vcf_entry.END}\t{vcf_entry.GENOTYPE}\t{vcf_entry.SVTYPE}\t{vcf_entry.SVLEN}\t{nSupport}\t{support_ids_notRef}', end="\n")

    def print_out_bed(vcf_entry):
        if vcf_entry.SVTYPE == "BND" or abs(vcf_entry.SVLEN) > args.minsize:
            start=min(vcf_entry.END, vcf_entry.POS)
            end=max(vcf_entry.END, vcf_entry.POS)
            if args.breakpoints_only and vcf_entry.SVTYPE != "BND" and vcf_entry.SVTYPE != "INS":
                print(f'{vcf_entry.CHROM}\t{start}\t{start}\t{vcf_entry.ID}\ts', end="\n")
                print(f'{vcf_entry.CHROM}\t{end}\t{end}\t{vcf_entry.ID}\te', end="\n")
            else:
                print(f'{vcf_entry.CHROM}\t{start}\t{end}\t{vcf_entry.ID}\tu', end="\n")

    def id2fam(family_table):
        id_2_fam = {}
        famfilehandl = open(family_table, "r")
        for l in famfilehandl:
            [smplID, fam] = l.rstrip("\n").split("\t")
            if smplID not in id_2_fam:
                id_2_fam[smplID] = fam
        famfilehandl.close()
        return(id_2_fam)

    def snf2_filtering(args, id_2_fam):
        if args.as_bed:
            pass
        else:
            print(f'#CHROM\tSTART\tEND\tGT\tSVTYPE\tSVLEN\tnSUPPORT\tSUPPORT(NO 0/0)', end="\n")
        for line in sys.stdin:
            # skip comments
            if line.startswith("#") and not line.startswith("##"):
                vcfHeader = vcf_header(line.rstrip("\n"))
            elif not line.startswith("#"):
                vcf_entry = vcf_line_pop(line.rstrip("\n"))
                nSupport = sum([int(x) for x in list(vcf_entry.SUPP_VEC)])
                if nSupport >= args.minsupp_calls:
                    support_ids = vcfHeader.get_sample_name(vcf_entry.SUPP_VEC, vcf_entry.all_GT, id_2_fam)
                    support_ids_notRef = [x for x in support_ids.split(", ") if ("0/0" not in x and "./." not in x)]
                    nSupport_ids_notRef = len(support_ids_notRef)
                    support_ids_notRef = ", ".join(support_ids_notRef)
                    if args.uniq_only:
                        if nSupport == 1:
                            if args.as_bed:
                                print_out_bed(vcf_entry)
                            else:
                                print_out(vcf_entry, nSupport, support_ids_notRef)
                    else:
                        if nSupport_ids_notRef != 0:
                            if args.as_bed:
                                print_out_bed(vcf_entry)
                            else:
                                print_out(vcf_entry, nSupport_ids_notRef, support_ids_notRef)

    def snf2_ont31(args, id_2_fam):
        #10kb
        _10kb = 10000
        print(f'#CHROM\tSTART\tEND\tSVTYPE\tSVLEN', end="\n")
        for line in sys.stdin:
            if line.startswith("#") and not line.startswith("##"):
                vcfHeader = vcf_header(line.rstrip("\n"))
            elif not line.startswith("#"):
                vcf_entry = vcf_line_pop(line.rstrip("\n"))
                support_ids_all = vcfHeader.get_sample_name(vcf_entry.SUPP_VEC, vcf_entry.all_GT, id_2_fam)
                # remove 0/0 GT
                support_ids = ", ".join([i for i in support_ids_all.split(", ") if "0/0" not in i])
                # filters: notFather, skipTRA, onlychrX, notOnlyinMother(not seen in proband)
                notFather = "F" not in support_ids
                skipTRA = vcf_entry.SVTYPE != "BND"
                onlychrXY = vcf_entry.CHROM == "X" or vcf_entry.CHROM == "Y"
                notOnlyMother = not("M" in support_ids and "P" not in support_ids) if args.family_table != "" else True
                sv_min_len_is_10kb = abs(vcf_entry.SVLEN) >= _10kb if vcf_entry.SVTYPE != "BND" else True
                if sv_min_len_is_10kb and notFather and skipTRA and onlychrXY and notOnlyMother and len(support_ids) > 0:
                    vcf_end_pos = vcf_entry.POS + vcf_entry.SVLEN if not isSV_translocation(vcf_entry.SVTYPE) else vcf_entry.TRA
                    print(f'{vcf_entry.CHROM}\t{vcf_entry.POS}\t{vcf_end_pos}\t{vcf_entry.SVTYPE}\t{vcf_entry.SVLEN}\t{support_ids}', end="\n")
            else:
                pass

    id_2_fam = id2fam(args.family_table) if args.family_table != "" else {}
    if args.ont31:
        if args.family_table == "":
            logger("Family relationship (--families) not provided", "warn")
        snf2_ont31(args, id_2_fam)
    else:
        snf2_filtering(args, id_2_fam)


def isSV_translocation(the_svtype):
    return (the_svtype == "BND" or the_svtype == "TRA")


def get_arguments():
    vcfu_help = """
    vcf_utils <command> [<args>]
            single      Extract genotype information for SV
                          none | --min-read-support | --filter-gt | --genotypes_only
            population  Parse information of Sniffles2 population merge result
                          none | --families | --min-support | --min-size | --uniq-only | --as-bed & --breakpoints-only | --ont-trios & --ont-trios-filter | --ont-31
            filterind   Filters SV from sniffles2 pop for a single indiviual (bcftools -s NAME pop.vcf)
                          none
            version    Version
    """
    # TODO: no command give no help
    parser = argparse.ArgumentParser(
             description="VCF utils",
             usage=vcfu_help
    )
    subparsers = parser.add_subparsers(help=vcfu_help, dest="command")

    # ############################################################################################ #

    # ############################################################################################ #
    # GenotypeSV
    genotype_sv_help = "Extracts the chr, position and genotype for SV"
    subparser_genotypesv = subparsers.add_parser("parsesv", help=genotype_sv_help)
    subparser_genotypesv.add_argument('-s', '--min-read-support', type=int, required=False, dest='minsupp_read', default=1,
                                   help='Min. read support for the SV calls, default = 1')
    subparser_genotypesv.add_argument('-f', '--filter-gt', type=str, required=False, dest='filer_gt', default="",
                                   help='Removed genotypes from output, for multiple,need to be comma separated.\nExample: -f 0/0,0/1, default = None')
    subparser_genotypesv.add_argument('-g', '--genotypes_only', action='store_true', required=False, dest='only_genotype',
                                   help='Only output the genotype counts, default = False')

    # ############################################################################################ #

    # ############################################################################################ #
    # Sniffles2 population
    population_help = "Perform analysis on Sniffles2 population-merges"
    subparser_population = subparsers.add_parser("population", help=population_help)
    subparser_population.add_argument('-f', '--families', type=str, required=False, dest='family_table', default="",
                                   help='File name or path of the file containing the relationship of a sampleID and its source. \
                                    Adds F for father, M for mother and P for proband in the results of the genotype')
    subparser_population.add_argument('-m', '--min-support', type=int, required=False, dest='minsupp_calls', default=1,
                                   help='Min. support for the SV calls (from SUPP_VEC), default = 1')
    subparser_population.add_argument('-s', '--min-size', type=int, required=False, dest='minsize', default=1,
                                   help='Min. absolute size of the event (excep for BDN), default = 1')
    subparser_population.add_argument('-u', '--uniq-only', action='store_true', required=False, dest='uniq_only',
                                   help='Show only those that appear in a single individual (from SUPP_VEC), default = False')
    subparser_population.add_argument('-b', '--as-bed', action='store_true', required=False, dest='as_bed',
                                   help='The output is bed format: chr start end (tab separated), default = False')
    subparser_population.add_argument('-k', '--breakpoints-only', action='store_true', required=False, dest='breakpoints_only',
                                   help='Write down only the breakpoints, e.g. start = A, end =A; so each SV is twice, once start, once end, default = False, Only works with "as-bed"')

    # Analysis for the paper
    subparser_population.add_argument('-1', '--ont-31', action='store_true', required=False, dest='ont31',
                                   help='Perform the analysis of the probands for the sniffles2 paper, default = False')

    # ############################################################################################ #

    # ############################################################################################ #
    # Version
    version_help = "Gives the version number"
    subparser_version = subparsers.add_parser("version", help=version_help)

    # ############################################################################################ #

    args = parser.parse_args()
    return args, vcfu_help # TODO: see above


def main():
    args, main_help = get_arguments()
    command = args.command

    # single
    if command == "parsesv":
        get_genotype_sv(args)
    # population
    elif command == "population":
        args.minsupp_calls = 1 if args.uniq_only else args.minsupp_calls # in uniq_only override minsupp_calls
        population_parse(args)
    # version/--version
    elif command == "version":
        the_version()
    # help
    else:
        print(main_help)

# main
if __name__ == '__main__':
    main()
