from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from Modules.SeqOrganizer import retTempDirectory

import sys, os, pysam
from collections import defaultdict

from subprocess import call, Popen, PIPE
from multiprocessing import cpu_count

from Modules.McGrathLargeVariant import ChimericCaller
from Modules.ElegansVarAnnTools import EVAT

class VariantCallerGenotyper():
    def __init__(self, analysis_mode, databaseID, discoverySamples, excludedSamples, genotypeSamples, population1, population2, delete = False):

        # Create dictionary for available analysis options. List contains alignment function, caller function, and then deletion function
        self.available_options = {}
        self.available_options['PatrickCaller'] = [self._patrickAligner, self._patrickCaller, None]

        # Save input values as attributes
        self.analysis_mode = analysis_mode
        self.databaseID = databaseID
        self.discoverySamples = set()
        self.excludedSamples = set()
        self.genotypeSamples = set()
        self.population1 = set()
        self.population2 = set()
        self.samples = {}

        # Perform error checking and populate good sample dictionary
        self._quality_control(discoverySamples, excludedSamples, genotypeSamples, population1, population2)

        self.delete = delete

    def identifyVariants(self):
        # Call alignment specific function for analysis mode:
        self.available_options[self.analysis_mode][0]()

        # Call variant caller specific function for analysis mode
        #self.available_options[self.analysis_mode][1]()
        
    def _quality_control(self, discoverySamples, excludedSamples, genotypeSamples, population1, population2):

        # Temporarary variables
        bad_flag = False
        badSamples = set()
        
        # Quality Control: Did the user give a valid analysis option?
        if not self._valid_option(self.analysis_mode):
            print(self.analysis_mode + ' is not a valid analysis option', file = sys.stderr)
            self._print_options()
            bad_flag = True
                
        # Quality Control: Is the database valid?
        if self.databaseID not in DI.retGenomeVersions():
            print('Reference database not found. Options are: ' + ','.join(DI.retGenomeVersions()), file = sys.stderr)
            bad_flag = True

        # Quality Control: Are the discoverySamples valid?
        for sample in discoverySamples:
            if sample not in SI.retSampleIDs():
                badSamples.add(sample)
            else:
                self.discoverySamples.add(sample)
                self.genotypeSamples.add(sample)
                self.samples[sample] = SI.retSampleData(sample)

        # Quality Control: Are the excludedSamples valid?
        if excludedSamples is not None:
            for sample in excludedSamples:
                if sample not in SI.retSampleIDs():
                    badSamples.add(sample)
                else:
                    self.excludedSamples.add(sample)
                    self.genotypeSamples.add(sample)
                    self.samples[sample] = SI.retSampleData(sample)

        # Quality Control: Are the genotypeSamples valid?
        if genotypeSamples is not None:
            for sample in genotypeSamples:
                if sample not in SI.retSampleIDs():
                    badSamples.add(sample)
                else:
                    self.genotypeSamples.add(sample)
                    self.samples[sample] = SI.retSampleData(sample)

        # population1 and population2 good
        if (population1 is None and population2 is not None) or (population1 is not None and population2 is None):
            print('If one population is provided the other population must be as well', file = sys.stderr)
            bad_flag = True
        if population1 is not None:
            for sample in population1:
                if sample not in SI.retSamplesIDs():
                    badSamples.add(sample)
                else:
                    self.genotypeSamples.add(sample)
                    self.samples[sample] = SI.retSampleData(sample)
                    self.population1.add(sample)

        if population2 is not None:
            for sample in population2:
                if sample not in SI.retSamplesIDs():
                    badSamples.add(sample)
                else:
                    self.genotypeSamples.add(sample)
                    self.samples[sample] = SI.retSampleData(sample)
                    self.population2.add(sample)
                            
        # Quality Control: If an error has occurred, print useful info and quit.
        if bad_flag or len(badSamples) > 0:
            if len(badSamples) > 0:
                print('The following bad samples were not found: ' + str(',').join(badSamples), file = sys.stderr)
                print('The available samples are: ' + str(',').join(SI.retSampleIDs()), file = sys.stderr)
            raise NotImplementedError
          
    def _patrickAligner(self):

        ref_files = DI.retDatabaseData(self.databaseID)
        ref_file = ref_files.refFile_z
        FNULL = open(os.devnull, 'w')
        sys.stdout = open(os.devnull, 'w')
        
        for sampleID in sorted(list(self.samples.keys())):
            
            sampleData = self.samples[sampleID]
            # If directory already exists, remove it
            if self.delete:
                FM.del_dir('PatrickDNA', self.databaseID, sampleID)

            # Create output file names
            all_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'all')
            unmapped_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'unmapped')
            discordant_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'discordant')
            inversion_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'inversion')
            duplication_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'duplication')
            clipped_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'clipped')
            chimeric_bm = FM.ret_asb('PatrickDNA', self.databaseID, sampleID, 'chimeric')
            
            bamfiles = []

            print('Aligning ' + sampleID, file = sys.stderr)
            print(all_bm, file = sys.stderr)
            if os.path.isfile(all_bm):
                print(sampleID + ' already exists. Skipping....', file = sys.stderr)
                continue
            for run in sampleData:

                # Create temporary file names
                tfile1 = FM.ret_temp_file()
                tfile2 = FM.ret_temp_file('.bam')
                tfile3 = FM.ret_temp_file('.bam')

                # Run BWA to map reads to reference
                if run.paired_flag:
                    if run.twofq:
                        print('   Processing paired fastq files ' + run.fqfile1.split('/')[-1] + ' and ' + run.fqfile2.split('/')[-1], file = sys.stderr)
                        bwa_command = ['bwa', 'mem', '-t', str(cpu_count()), '-R', run.RG_info, '-M', ref_file, run.fqfile1, run.fqfile2]
                        #print('   ' + ' '.join(bwa_command), file = sys.stderr)
                        call(bwa_command, stdout = open(tfile1, 'w'), stderr = FNULL)
                    else:
                        print('   Processing interleaved fastq file ' + run.fqfile1.split('/')[-1], file = sys.stderr)
                        bwa_command = ['bwa', 'mem', '-t', str(cpu_count()), '-p', '-R', run.RG_info,'-M', ref_file, run.fqfile1]
                        #print('   ' + ' '.join(bwa_command))
                        call(bwa_command, stdout = open(tfile1, 'w'), stderr = FNULL)
                else:
                    print('   Processing SE fastq file ' + run.fqfile1.split('/')[-1], file = sys.stderr)
                    bwa_command = ['bwa', 'mem', '-t', str(cpu_count()), '-R', run.RG_info, '-M', ref_file, run.fqfile1]
                    #print('   ' + ' '.join(bwa_command), file = sys.stderr)
                    call(bwa_command, stdout = open(tfile1, 'w'), stderr = FNULL)

                # Sort sam file and convert to bam file
                print('   Sorting file...', file = sys.stderr)
                p1 = Popen(['samtools', 'view', '-bh', '-@', str(cpu_count()), tfile1], stdout=PIPE)
                p2 = Popen(['samtools', 'sort','-o', tfile2, '-@', str(cpu_count()), '-'], stdin = p1.stdout, stderr = FNULL)
                p2.communicate()

                # Remove duplicates
                print('   Removing duplicates...', file = sys.stderr)
                call(['samtools','rmdup','-s', tfile2, tfile3], stderr = FNULL)
                #call(['java', '-jar','/usr/local/share/java/picard.jar','SortSam', 'I=' + tfile1, 'O=' + tfile2, 'SORT_ORDER=coordinate'])

                # Remove temporary files and append last file to array to merge
                call(['rm', '-f', tfile1, tfile2, 'temp.txt'])
                bamfiles.append(tfile3)

            # Merge bamfiles together
            print('Merging bamfiles together and indexing master file...', file = sys.stderr)
            if len(bamfiles) > 1:
                call(['samtools', 'merge', '-f', all_bm] + bamfiles)
            else:
                call(['mv', bamfiles[0], all_bm])
            # Index masterfile
            call(['samtools', 'index', all_bm])

            # Remove remaining temporary files
            call(['rm', '-f'] + bamfiles)

            # Open master bam file
            align_file = pysam.AlignmentFile(all_bm) 
            # Open bam files that will hold reads split by alignment type
            unmapped = pysam.AlignmentFile(unmapped_bm, mode = 'wb', template = align_file) 
            discordant = pysam.AlignmentFile(discordant_bm, mode = 'wb', template = align_file)
            inversion = pysam.AlignmentFile(inversion_bm, mode = 'wb', template = align_file)
            duplication = pysam.AlignmentFile(duplication_bm, mode = 'wb', template = align_file)
            clipped = pysam.AlignmentFile(clipped_bm, mode = 'wb', template = align_file)
            chimeric = pysam.AlignmentFile(chimeric_bm, mode = 'wb', template = align_file)

            # Go through all reads and process them into appropriate categories
            print('Splitting reads based upon their alignment...', file = sys.stderr)
            for read in align_file.fetch(until_eof=True):
                if read.is_paired:
                    # Both reads are unmapped
                    if read.is_unmapped and read.mate_is_unmapped:
                        unmapped.write(read)
                    # One read is unmapped
                    elif read.is_unmapped or read.mate_is_unmapped:
                        discordant.write(read)                
                    # Chromosome fusion
                    elif read.reference_id!=read.next_reference_id:
                        discordant.write(read)
                    # Inversion
                    elif read.is_reverse == read.mate_is_reverse:
                        inversion.write(read)
                    # Duplication
                    elif ((read.pos < read.mpos and read.is_reverse) or (read.pos > read.mpos and read.mate_is_reverse)) and abs(read.isize) > 102:
                        duplication.write(read)
                else:
                    if read.is_unmapped:
                        unmapped.write(read)
                # Clipped
                if read.cigarstring is not None:
                    for pair in read.cigartuples:
                        if pair[0] == 4 and pair[1] > 4:
                            clipped.write(read)
                            break
                        elif pair[0] == 5 and pair[1] > 4:
                            clipped.write(read)
                            break
                # Chimeric
                if read.has_tag('SA'):
                    chimeric.write(read)

            align_file.close()
            unmapped.close()
            discordant.close()
            inversion.close()
            duplication.close()
            clipped.close()
            chimeric.close()

            # Index new bamfiles
            call(['samtools','index', unmapped_bm])
            call(['samtools','index', discordant_bm])
            call(['samtools','index', inversion_bm])
            call(['samtools','index', duplication_bm])
            call(['samtools','index', clipped_bm])
            call(['samtools','index', chimeric_bm])
 
            qualimap_dir = FM.ret_qualimap_dir('PatrickDNA', self.databaseID, sampleID)
            # Hack to deal with inability of qualimap to deal with spaces in file name
            print('\t'.join(['/Users/pmcgrath7/qualimap_v2.2.1/qualimap', 'bamqc', '-bam', all_bm.replace(os.getcwd() + '//', ''), '-c', '-outdir', qualimap_dir.replace(os.getcwd() + '//', '')]), file = sys.stderr)
            call(['/Users/pmcgrath7/qualimap_v2.2.1/qualimap', 'bamqc', '-bam', all_bm.replace(os.getcwd() + '//', ''), '-c', '-outdir', qualimap_dir.replace(os.getcwd() + '//', '')])

    def _patrickCaller(self):
        ref_files = DI.ret_database_data(self.genomeVersion, self.species)
        ref_file = ref_files.ref_file
        unzipped_ref_file = ref_files.unzipped_ref_file

        # Identify small variants
        
        samples = []
        samples_bm = []
        for sampleID in sorted(list(self.good_samples.keys())):
            print('Identifying small variants for ' + sampleID, file = sys.stderr)
            projectID = self.good_samples[sampleID][0].projectID
            samples.append(sampleID)
            samples_bm.append(FM.ret_asb('PatrickDNA', self.species, self.genomeVersion, projectID, sampleID, 'all'))
        sample_vcf = FM.ret_vcf_file('Polys:SNVs', self.species, self.genomeVersion, samples, 'None', samples)
        print('Identifying small variants from the following samples: ' + ', '.join(samples), file = sys.stderr)
        #call(['freebayes', '-f', unzipped_ref_file, '-C', '6', '-F', '0.8'] + samples_bm, stdout = open(sample_vcf, 'w'))
        with open(sample_vcf) as f, open(sample_vcf.replace('.vcf', '.filtered.vcf'), 'w') as g:
            for line in f:
                if line[0] == '#':
                    print(line, file = g)
                    if line[0:6] == '#CHROM':
                        ref_index = line.rstrip().split('\t').index(self.reference_sample)
                else:
                    if line.split('\t')[ref_index][0:3] == '0/0':
                        print(line.rstrip(), file = g)

        call(['snpeff', 'eff', self.genomeVersion, sample_vcf.split('ElegansIndelSequencer/')[1].replace('.vcf', '.filtered.vcf')], stdout = open(sample_vcf.replace('.vcf', '.filtered.annotated.vcf'), 'w'))
        
        #call(['bgzip', sample_vcf])
        #call(['tabix', '-p', 'vcf', sample_vcf + '.gz'])

        refprojectID = self.good_samples[self.reference_sample][0].projectID
        all_bams = {}
        chimeric_bams = {}
        for sampleID in sorted(list(self.good_samples.keys())):
            projectID = self.good_samples[sampleID][0].projectID
            all_bams[sampleID] = FM.ret_asb('PatrickDNA', self.species, self.genomeVersion, projectID, sampleID, 'all')
            chimeric_bams[sampleID] = FM.ret_asb('PatrickDNA', self.species, self.genomeVersion, projectID, sampleID, 'chimeric')
            #Create chimeric caller object and identify large deletions/variants
        sample_vcf2 = FM.ret_vcf_file('Polys:LargeVar', self.species, self.genomeVersion, samples, 'None', samples)
        cc = ChimericCaller(chimeric_bams, all_bams, list(self.good_samples.keys()), list(self.good_samples.keys()), self.reference_sample, ref_file, sample_vcf2)
        cc.identifyChimericVariants()
        call(['snpeff', 'eff', self.genomeVersion, sample_vcf2.split('ElegansIndelSequencer/')[1]], stdout = open(sample_vcf2.replace('.vcf', '.annotated.vcf'), 'w'))
        
        summarized_data = FM.ret_vcf_file('Polys:Summarized', self.species, self.genomeVersion, samples, 'None', samples).replace('.vcf', '.txt')
        EVAT([sample_vcf.replace('.vcf', '.filtered.annotated.vcf'),sample_vcf2.replace('.vcf', '.annotated.vcf')], summarized_data, self.population1, self.population2)
           

    def _valid_option(self, option):
        if option in self.available_options:
            return True
        else:
            return False
        
    def _print_options(self):
        print('Available analysis options are: ' + ', '.join(list(self.available_options.keys())), file = sys.stderr)
        print('', file = sys.stderr)
        

