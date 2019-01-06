########################################################################
#
#   Module DATAUNITE.
#       To unite all tables regarding long non-coding RNAs data in Mango.
#       Current data:
#       o   Structure stability.
#       o   Orthologues predictions 40 qcov_hsp.
#       o   Orthologues predictions 50 qcov_hsp.
#       o   Differential expressed lnc-RNAs.
#       o   Differential expressed total RNAs.
#       o   FPKMs of lnc-RNAs
#       o   FPKMs of total RNAs
#       
#   
#   Author:     Javier Montalvo-Arredondo.
#   Contact:    Biotechnology divison at
#               Universidad Autonoma Agraria Antonio Narro,
#               Calzada Antonio Narro S/N,
#               San Buenavista, Saltillo, Coahuila, Mexico.
#               buitrejma[at]gmail.com
#
#       TIDY TABLE STRUCTURE.
#       o   Factor: [RNA_with_ORF, RNA_whitout_ORF]
#       o   Factor: [has_orth_50, has_orth_40, has_no_orth
#       o   Factor: [stable_struct, unstable_struct]
#                   #!#!#!# [WARNING!] - PERFORM RNAfold for every RNA; select those structures stables and unstables
#                           based on predicted CDS sequence stability.
#       o   Factor: [more_stable_struct_than_pred_CDS, equal_stable_struct_then_pred_CDS, less_stable_struct_than_pred_CDS, is_a_pred_CDS]
#       o   Factor: [is_DEG, isnt_DEG]
#       o   Variable: [FPKMs]
#       o   Variable: [struct_DGEs] # [HERE I AM!]
#       o   Variable: [norm_struct_DGEs]
#       o   Variable: [struct_DGEs_pval]
#       o   Variable: [percent_identity]
#       o   Variable: [percent_alignment]
#       o   Variable: [qcov_hsp]
#       o   Variable: [au_gc_ratio]
#
########################################################################
#
#   
#

#   [TESTING] - [OK]
def parsing_args():
    """
        Returns a hash table where the key is a tag and the
        value is the argument.
        for example:
        the following are passed arguments as input:
            -i infile -o outfile -tmp temp_file
        the following is the expected hash table:
            hash =  {
                    "-i":"infile",
                    "-o":"outfile",
                    "-tmp":"tmp_file"
                    }
    """
    from sys import argv, exit
    args_and_tags = argv[1:]
    parsed_args = {}
    if len(args_and_tags) % 2 != 0:
        print "[FATAL ERROR!] - Bad arguments number"
        exit(1)
    else:
        for i in xrange(len(args_and_tags)):
            if i % 2 == 0:
                tag = args_and_tags[i]
                if tag[0] != "-":
                    print "[FATAL ERROR!] - Tag does not begin with \"-\""
                    exit(1)
            else:
                 parsed_args[tag] = args_and_tags[i].split(" ")\
                                    if len(args_and_tags[i].split(" ")) > 1\
                                    else args_and_tags[i]

    return parsed_args

def IF_serialise(line, values_idx):
    """
        Takes the array values_idx an it will serialise the data
        in line which data index is equal to the indexes of the array
        values_idx.
    """
    array_out = []
    for idx in values_idx:
        array_out.append(line[idx])

    return array_out

def IF_hashing(filename, key_idx, values_idx):
    """
        Transform a plain-text data table into a hash table.
        The key and value or values are defined in parsed args.
        for instance:
            args[--rnas-with-orfs] stores the data table.
            args[--rnas-with-orfs-key] stores the number of the column
            of the data table that will be used as key in the hash table.
            args[--rnas-with-orfs-values] stores the number o numbers
            of the columns of the data table from which data is obtained.
    """
    FH = open(filename, "r")
    hash_table = {}
    array = []
    key_idx = int(key_idx)
    values_idx = [int(x)-1 for x in values_idx] # [BUG1] - Values indexes were not 0-based indexes.
    len_values_idx = len(values_idx)
    key_idx -= 1
    for line in FH:
        line = line.strip("\n")
        line = line.split("\t")
        if values_idx == "None":
            hash_table[line[key_idx]] = None
        else:
            hash_table[line[key_idx]] = IF_serialise(line, values_idx)

    hash_table["len_values_idx"] = len_values_idx
    FH.close()
    
    return hash_table

def IF_arraying(filename):
    """
        Extract lines from filename file and stores them in an array.
    """
    FH = open(filename, "r")
    array = []
    for line in FH:
        line = line.strip("\n")
        array.append(line)

    return array

def factorise(parsed_args):
    """
        This function loads the data tables of lncrnas project
        and extracts specific data for variables and factors
        for every gene ID.
        Every factor and variable produced is a hash table and is
        stored into another hash table where key is the factor name or
        variable name.
        It requires the parsed arguments produced by parsing_args().
    """
    from sys import stdout as out
    # LOADING FACTOR1: [RNA_with_ORF, RNA_without_ORF]
    args = parsed_args
    #   [WARNING!] - IF_function means that this function
    #   is for internal purposes.
    array_all_headers = IF_arraying(args["--all-headers"])
    array_ORFS_headers = IF_arraying(args["--rnas-with-orfs"])
    array_NOORFS_headers = IF_arraying(args["--rnas-without-orfs"])
    array_has_orth_50 = IF_arraying(args["--has-orth-50"])
    array_has_orth_40 = IF_arraying(args["--has-orth-40"])
    array_stable_struct = IF_arraying(args["--stable-struct"])
    array_DEG_headers = IF_arraying(args["--is-DEG"])
    hash_exp = IF_hashing(args["--exp"],\
                        args["--exp-key"],\
                        args["--exp-val"])
    hash_stat_dgs = IF_hashing(args["--stats-dgs"],\
                            args["--stats-dgs-key"],\
                            args["--stats-dgs-val"])
   # array_has_noorth = IF_arraying(args["--has-no-orth"])
    f1 = {}
    total = float(len(array_all_headers))
    counter = 0.0
    for header in array_all_headers:
        #   [PROCEDURE!] - Selecting RNA with ORFS & RNA without ORFs from total RNAs.
        if header in array_ORFS_headers:
            f1[header] = ["RNA_with_ORF"]
        elif header in array_NOORFS_headers:
            f1[header] = ["RNA_without_ORF"]
        else:
            f1[header] = ["NULL"]
        #   [PROCEDURE] - Selecting RNAs with orthologues from RNAs without ORFs.
        if header in array_has_orth_50 and header in array_NOORFS_headers:
            f1[header].append("has_orth_50")
        elif header in array_has_orth_40 and header in array_NOORFS_headers:
            f1[header].append("has_orth_40")
        elif header in array_NOORFS_headers:
            f1[header].append("has_no_orth")
        else:
            f1[header].append("NULL")
        #   [PROCEDURE] - Selecting RNAs without ORFs with stable structure.
        if header in array_NOORFS_headers:
            if header in array_stable_struct:
                f1[header].append("stable_struct")
            else:
                f1[header].append("unstable_struct")
        else:
            f1[header].append("NULL")
        #   [PROCEDURE] - Selecting RNAs without ORFs with a DEG expression.
        if header in array_DEG_headers:
            f1[header].append("is_DEG")
        else:
            f1[header].append("is_not_DEG")
        #   [FATAL ERROR!] - It did not work - [FATAL ERROR!] # [BUG1] - values_idx were not 0-based index
        if header in hash_exp.keys():
            f1[header].extend(hash_exp[header])
        else:
            f1[header].extend(["NULL" for x in xrange(hash_exp["len_values_idx"])])

        if header in hash_stat_dgs.keys():
            f1[header].extend(hash_stat_dgs[header])
        else:
            f1[header].extend(["NULL" for x in xrange(hash_stat_dgs["len_values_idx"])])
        
        out.write("Progress precentage: %d%% \r" % ((counter / total) * 100))
        out.flush()
        counter +=1
    
    out.write("Progress percentage: 100%\n[SUCCESS!] - Task finished")
    out.write("\n")
    return f1


def write_out(f1, output):
    """
        Writes out f1 to file: f1 is a hash table with all datat united.
    """
    FHOUT = open(output, "w+")
    for key, value in f1.iteritems():
        str_out = key\
                + "\t"\
                + "\t".join(value)\
                + "\n"
        FHOUT.write(str_out)

    FHOUT.close()

