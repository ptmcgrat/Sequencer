import argparse, sys

import Modules.DatabaseWorker as DBW
import Modules.ReadsWorker as RW
import Modules.AnalysisWorker as AnW
import Modules.VariantCallerWorker as VCW

from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI



parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='Available Commands', dest='command')

db_parser = subparsers.add_parser('Databases', help='Download or remove databases from WormBase or manipulate local databases')
db_parser.add_argument('-d', '--download', type = str, nargs = '+', help = 'List of databases (<species>,<wormbase>) to download from WormBase. If no option is provided, the list of available databases will be printed out')
db_parser.add_argument('-a', '--add_local_databases', type = str, help = 'An excel files containing info on local database you would like to add to the database.')
db_parser.add_argument('-r', '--remove', type = str, nargs = '+', help = 'Remove databases permanently.')

reads_parser = subparsers.add_parser('Reads', help='Download reads from SRA/ENA or manipulate newly generated reads')
reads_parser.add_argument('-d', '--download', type = str, nargs = '+', help = 'List of projects you would like to download from the SRA/ENA')
reads_parser.add_argument('-a', '--add_local_reads', type = str, help = 'List of excel files containing info on local reads you would like to add to the database.')
reads_parser.add_argument('-r', '--remove', type = str, nargs = '+', help = 'Remove data permanently. Can be used with the -s option')

VariantCaller_parser = subparsers.add_parser('VariantCaller', help = 'This command calls potential genetic variants for the samples you input. It performs no quality control or genotyping with respect to multiple samples')
VariantCaller_parser.add_argument('AnalysisType', type = str, help = 'Analysis you would like to run.')
VariantCaller_parser.add_argument('DatabaseID', type = str, help = 'Version of genome you want to use. Assumes appropriate files are installed')
VariantCaller_parser.add_argument('DiscoverySamples', type = str, nargs = '+', help = 'List of sample ids you would like to discover large variants from')
VariantCaller_parser.add_argument('-es', '--ExcludedSamples', type = str, nargs = '+', help = 'List of sample ids you would like to exclude large variants from')
VariantCaller_parser.add_argument('-gs', '--GenotypeSamples', type = str, nargs = '+', help = 'List of sample ids you would like to genotype large variants from')
VariantCaller_parser.add_argument('-p1', '--PopulationOne', type = str, nargs = '+', help = 'List of sample ids that belong to population 1')
VariantCaller_parser.add_argument('-p2', '--PopulationTwo', type = str, nargs = '+', help = 'List of sample ids that belong to population 2')
VariantCaller_parser.add_argument('-d', '--Delete', action = 'store_true', help = 'Delete existing alignment data before running')



AnalyzeVCF_parser = subparsers.add_parser('AnalyzeVCF', help='Align sequencing reads to a given database')
AnalyzeVCF_parser.add_argument('AnalysisType', type = str, help = 'Analysis you would like to run. Available options are: ')
AnalyzeVCF_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
AnalyzeVCF_parser.add_argument('IDs', type = str, nargs = '+',  help = 'List of projects or strains you would like to be aligned. <all> = all strains')
AnalyzeVCF_parser.add_argument('-s', '--species', help = 'Specify which species to align to')

AnalyzeVCF_parser.add_argument('-p1', '--PopulationOne', type = str, nargs = '+', help = 'List of sample ids that belong to population 1')
AnalyzeVCF_parser.add_argument('-p2', '--PopulationTwo', type = str, nargs = '+', help = 'List of sample ids that belong to population 2')


args = parser.parse_args()

if args.command is None:
    parser.print_help()

if args.command == 'Databases':
    try:
        db_obj = DBW.DBMaker(args.download, args.add_local_databases, args.remove)
    except NotImplementedError:
        print('Incorrect use of Databases script', file = sys.stderr)
        db_parser.print_help()

if args.command == 'Reads':
    try:
        reads_obj = RW.ReadWorker(args.download, args.add_local_reads, args.remove)
    except NotImplementedError:
        print('Incorrect use of Reads script', file = sys.stderr)
        reads_parser.print_help()        
    

if args.command == 'VariantCaller':

    try:
        caller_obj = VCW.VariantCallerGenotyper(args.AnalysisType, args.DatabaseID, args.DiscoverySamples,  args.ExcludedSamples, args.GenotypeSamples, args.PopulationOne, args.PopulationTwo, delete = args.Delete)
    except NotImplementedError:
        print('Error in input values to program. Quitting', file = sys.stderr)
        #VariantCaller_parser.print_help()
        sys.exit(1)

    caller_obj.identifyVariants()
    
if args.command == 'AnalyzeProject':
    print('Aligning Reads...', file = sys.stderr)
    if args.AnalysisType not in analysis:
        print('Not a good analysis option. Options are: ' + str(analysis), file = sys.stderr)
    aln_obj = AnW.AnalysisMaker(args.AnalysisType, args.Genome_version, args.IDs, args.species)
    


if args.command == 'VariantAnalyzer':

    if args.large:
        pass
    else:
        small_caller_obj = VCW.SNVCallerGenotyper(args.Species, args.Genome_version)
        small_caller_obj.analyzeVariants(args.IncludedSamples, args.ExcludedSamples, args.AnalyzedSamples)
    if args.small:
        pass
    else:
        large_caller_obj = VCW.LargeVarCallerGenotyper(args.Species, args.Genome_version)
        large_caller_obj.analyzeVariants(args.IncludedSamples, args.ExcludedSamples, args.AnalyzedSamples)
