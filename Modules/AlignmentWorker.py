from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from subprocess import call, Popen, PIPE
from multiprocessing import cpu_count
import sys, os, pysam


class AlignDNAMaker():
    def __init__(self, alignment_mode, genome_version, IDs, species, unlink):
        self.alignment_mode = alignment_mode
        self.genome_version = genome_version
        self.species = species
        self.IDs = IDs
        self.unlink = unlink

        bad_flag = False
        inc_species = set()
        
        for ID in self.IDs:
            if not SI.valid_project(ID) and not SI.valid_sample(ID):
                print(ID + ' is not a valid projectID or a project sampleID. Quitting...', file = sys.stderr)
                bad_flag = True
            elif SI.valid_project(ID):
                for sample in SI.ret_samples(ID):
                    for fq in SI.ret_sample_data(sample, ID):
                        inc_species.add(fq.species)

            elif SI.valid_sample(ID):
                for fq in SI.ret_sample_data(sample):
                    inc_species.add(fq.species)

        if species is None and len(inc_species) > 2:
            print('Samples have multiple species associated with them: ' + str(inc_species) + '. Quitting...', file = sys.stderr)
            bad_flag = True

        if species is None:
            self.species = next(iter(inc_species))

        if not DI.valid_database(self.genome_version, self.species):
            print('Reference database not found. Species: ' + self.species + '\tGenomeVersion: ' + self.genome_version, file = sys.stderr)
            bad_Flag = True
            
        if bad_flag:
            sys.exit()

        self.ref_files = DI.ret_database_data(self.genome_version, self.species)

    def align(self):
        #if self.split:
        #    print('Split mode chosen, will not realign fastq files to create bam files', file = sys.stderr)
        #    return

        if self.unlink:
            print('Unlink mode chosen, will unlink dropbox folder from local drive', file = sys.stderr)
            print('WARING: UNLINK MODE NOT IMPLEMENTED', file = sys.stderr)
            return

        for ID in self.IDs:
            if SI.valid_project(ID):
                samples = SI.ret_samples(ID)
                projectID = ID
                
            elif SI.valid_sample(ID):
                samples = [ID]
                SI.ret_sample_data(sample)
                projectID = SI.ret_sample_data(sample)[0].projectID
                
            if self.alignment_mode == 'PatrickRNA':
                self._alignPatrickRNA(projectID, samples)
            elif self.alignment_mode == 'PatrickDNA':
                self._alignPatrickDNA(projectID, samples)

    def _alignPatrickDNA(self, projectID, samples):

 

    def alignPatrickRNA(self, projectID, sampleID, fqs):

        FM.del_dir('PatrickRNA', self.species, self.genome_version, projectID, sampleID)
        bamfiles = []
        ref_file = self.ref_files.ref_file
        for fq in fqs:
            print(fq.fqfile1)
            tfile1 = FM.ret_temp_file()
            print('Aligning with bwa to ' + tfile1, file = sys.stderr)

            if fq.paired_flag:
                if fq.two_fq_flag:
                    print(' '.join(['bwa', 'mem', '-t', str(cpu_count()), '-R', fq.read_group, '-M', ref_file, fq.fqfile1, fq.fqfile2]), file = sys.stderr)
                    call(['bwa', 'mem', '-t', str(cpu_count()), '-R', fq.read_group, '-M', ref_file, fq.fqfile1, fq.fqfile2], stdout = open(tfile1, 'w'))
                else:
                    print(' '.join(['bwa', 'mem', '-t', str(cpu_count()), '-p', '-R', fq.read_group, '-M', ref_file, fq.fqfile1]))
                    call(['bwa', 'mem', '-t', str(cpu_count()), '-p', '-R', fq.read_group,'-M', ref_file, fq.fqfile1], stdout = open(tfile1, 'w'))
                    #call('bwa mem -t ' + str(cpu_count()) + ' -p -R ' + fq.RG_info + ' -M Databases//c_elegans/WS220/c_elegans.WS220.genomic.fa.gz Reads/reads.fq', stdout = open(tfile1, 'w'), shell = True)

            else:
                print(' '.join(['bwa', 'mem', '-t', str(cpu_count()), '-R', fq.read_group, '-M', ref_file, fq.fqfile1]), file = sys.stderr)
                call(['bwa', 'mem', '-t', str(cpu_count()), '-R', fq.read_group, '-M', ref_file, fq.fqfile1], stdout = open(tfile1, 'w'), stderr = open('bam_mem.txt', 'w'))
            tfile2 = FM.ret_temp_file('.bam')
            print('Sorting with samtools to ' + tfile2, file = sys.stderr)
            p1 = Popen(['/usr/local/bin/samtools', 'view', '-bh', tfile1], stdout=PIPE)
            p2 = Popen(['/usr/local/bin/samtools', 'sort','-o', tfile2, '-@', str(cpu_count()), '-'], stdin = p1.stdout)
            p2.communicate()
            call(['rm', '-f', tfile1, 'temp.txt', 'bam_mem.txt'])
            bamfiles.append(tfile2)
            
        out_file = FM.ret_asb('PatrickRNA', self.species, self.genome_version, projectID, sampleID, 'all')
        
        if len(bamfiles) > 1:
            call(['/usr/local/bin/samtools', 'merge', '-f', out_file] + bamfiles)
        else:
            print(['mv', bamfiles[0], out_file])
            call(['mv', bamfiles[0], out_file])
        call(['/usr/local/bin/samtools', 'index', out_file])
        call(['rm', '-f'] + bamfiles)
                    
    def merge_groups(self):
        if self.group is False:
            return
        [fg,sg] = self.group.split('_')
        trues = []
        falses = []
        bf_True, bf_False = FM.ret_asb(self.group, self.genome_version, 'discordant')
        for sample, grouping in SI.groupings[self.group].items:
            if grouping == True:
                trues.append(FM.ret_asb(sample, self.genome_version, 'discordant'))
            elif grouping == False:
                falses.append(FM.ret_asb(sample, self.genome_version, 'discordant'))
        call(['samtools', 'merge', bf_True, trues])
        call(['samtools', 'merge', bf_False, falses])

    
