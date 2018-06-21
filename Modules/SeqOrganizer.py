#Changed
import os, inspect, sys, pandas, shutil
from subprocess import call, Popen, PIPE

import os, sys, random, inspect, xlrd, pandas, shutil, multiprocessing
from collections import defaultdict
import pdb

#Deterimine absolute path of where this program is installed
installed_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/'

#Read in file containing necessary locations of files
definitions = {}
with open(installed_directory + 'Definitions.txt') as f:
    for line in f:
        key = line.split('=')[0].replace(' ','')
        data = line.rstrip().split('=')[1].replace(' ','').replace('~', os.path.expanduser('~'))
        definitions[key] = data
        
def retTempDirectory():
    return definitions['TempDirectory']

# Where data will be stored
output_directory = installed_directory + '/Output/'
polymorphism_directory = installed_directory + '/Polymorphisms/'
analysis_directory = installed_directory + '/Analysis/'

class RefObj():
    def __init__(self, databaseID, species, basedir, gtf_flag = False):
        self.databaseID = databaseID
        self.species = species
        self.basedir = basedir
        if gtf_flag == 'True':
            self.gtf_flag = True
        elif gtf_flag == 'False':
            self.gtf_flag = False
        else:
            self.gtf_flag = gtf_flag
        self.refFile_z = basedir + 'Genome.fa.gz'
        self.refFile_uz = basedir + 'Genome.fa'
        if gtf_flag:
            self.annotationFile = basedir + 'genes.gtf'
        else:
            self.annotationFile = basedir + 'genes.gff'
        self.transcriptFile = basedir + 'Transcript.fa.gz'
        self.kallistoFile = self.transcriptFile.split('fa.gz')[0] + 'idx'

    def retDataRow(self):
        out = {}
        out['DatabaseID'] = self.databaseID
        out['Species'] = self.species
        out['BaseDir'] = self.basedir
        out['GTFFlag'] = self.gtf_flag
        return out
    
class RefDatabase():
    def __init__(self):
        self.dbDirectory = installed_directory + '/Databases/'
        if not os.path.exists(self.dbDirectory):
            os.makedirs(dbDirectory)
        self.refDatabase =  self.dbDirectory + '/AvailableDatabases.xlsx'
        self.dt = pandas.read_excel(self.refDatabase)

    def addDatabase(self, newDatabase):
        # This file takes in a dictionary containing files and info necessary for a new database
        
        # 1. Reread data file to make sure it hasn't changed
        self.dt = pandas.read_excel(self.refDatabase)

        needed_data = ['DatabaseID', 'Species', 'GenomicFile', 'AnnotationFile', 'TranscriptFile']

        # 2. Check to make sure data is included
        for datum in needed_data:
            if datum not in newDatabase:
                print('Data dictionary must include ' + datum + ' key to add database. Exiting....', file = sys.stderr)
                raise KeyError

        # 3. Check if DatabaseID is already present
        if len(self.dt[self.dt.DatabaseID == newDatabase['DatabaseID']]) != 0:
            print(newDatabase['DatabaseID'] + ' already present in the database. Cannot add database with the same DatabaseID', file = sys.stderr)
            raise KeyError

        # 4. baseDir where data will be stored
        baseDir = self.retDatabaseDir(newDatabase['DatabaseID'], newDatabase['Species'])

        # 5. Determine if gtf or gff file is present
        if '.gtf' in newDatabase['AnnotationFile']:
            gtf_flag = True
        elif '.gff' in newDatabase['AnnotationFile']:
            gtf_flag = False
        else:
            print('Not sure what format ' + newDatabase['AnnotationFile'] + ' is. Should end with .gff or .gtf', file = sys.stderr)
        
        tobj = RefObj(newDatabase['DatabaseID'], newDatabase['Species'], baseDir, gtf_flag)
        
        # 6. Move files to appropriate directories
        if not os.path.isfile(newDatabase['GenomicFile']):
            print(newDatabase['GenomicFile'] + ' not found')
            raise FileNotFoundError
        if '.gz' not in newDatabase['GenomicFile']:
            call(['gzip', '-c', newDatabase['GenomicFile']], stdout = open(tobj.refFile_z, 'w'))
        else:
            shutil.copyfile(newDatabase['GenomicFile'], tobj.refFile_z)

        if not os.path.isfile(newDatabase['AnnotationFile']):
            print(newDatabase['AnnotationFile'] + ' not found')
            raise FileNotFoundError
        if '.gz' in newDatabase['AnnotationFile']:
            call(['gunzip', newDatabase['AnnotationFile']])
            newDatabase['AnnotationFile'] = newDatabase['AnnotationFile'].replace('.gz','')
 
        shutil.copyfile(newDatabase['AnnotationFile'], tobj.annotationFile)

        if not os.path.isfile(newDatabase['TranscriptFile']):
            print(newDatabase['TranscriptFile'] + ' not found')
            raise FileNotFoundError
        shutil.copyfile(newDatabase['TranscriptFile'], tobj.transcriptFile)

        if 'AnalysisDirectory' in newDatabase:
            shutil.move(newDatabase['AnalysisDirectory'], baseDir)
  
        
        # 7. Process files
        print('Indexing and processing reference files', file = sys.stderr)
        call(['gunzip', tobj.refFile_z])
        call(['bgzip', '-@', '4', '-c', tobj.refFile_uz], stdout = open(tobj.refFile_z, 'w'))

        call(['bwa', 'index', tobj.refFile_z])

        call(['samtools', 'faidx', tobj.refFile_z])
        call(['samtools', 'faidx', tobj.refFile_uz])

        call(['kallisto', 'index', tobj.transcriptFile, '-i', tobj.kallistoFile])

        # Add to snpeff file
        snpeffDirectory = definitions['SnpEffDataDir'] + '/' + tobj.databaseID + '/'
        try:
            os.makedirs(snpeffDirectory)
        except FileExistsError:
            pass
        shutil.copyfile(tobj.refFile_uz, snpeffDirectory + 'sequences.fa')
        with open(definitions['SnpEffConfigFile']) as f:
            contents = f.readlines()
        with open(definitions['SnpEffConfigFile'], 'w') as f:
            i = contents.index('# Non-standard Databases\n')
            output = contents[:i+3] + ['# ' + tobj.species + ' ' + tobj.databaseID + '\n', tobj.databaseID + '.genome : ' + tobj.databaseID + '\n', '\n'] + contents[i+3:]
            f.write(''.join(output))

        if tobj.gtf_flag:
            shutil.copyfile(tobj.annotationFile, snpeffDirectory + 'genes.gtf')
            call(['snpeff', 'build', '-gtf22', tobj.databaseID])

        else:
            shutil.copyfile(tobj.annotationFile, snpeffDirectory + 'genes.gff')
            call(['snpeff', 'build', '-gff3', tobj.databaseID])

        """
        cur_path = os.path.realpath(__file__)
        script_path = '/'.join(cur_path.split('/')[:-2])+'/Scripts'
        file_prefix = data['Species']+'.PRJNA13758'+'.'+data['Version']
        gtf_file = data['GtfFile']
        ss_file = file_prefix+'.ss'
        exon_file = file_prefix+'.exon'
        p1 = Popen(['gzip','-cd', data_dir+gtf_file], stdout=PIPE)
        p2 = Popen([script_path+'/hisat2_extract_splice_sites.py', '-'], stdin = p1.stdout,stdout = open(data_dir +ss_file, 'w'))
        p2.communicate()
        p1 = Popen(['gzip','-cd', data_dir+gtf_file], stdout=PIPE)
        p2 = Popen([script_path+'/hisat2_extract_exons.py', '-'], stdin = p1.stdout,stdout = open(data_dir +exon_file, 'w'))
        p2.communicate()
        l_ref = data_dir +data['GenomicFile']
        l_hisat2_index = data_dir + file_prefix + '.hisat2_index'
        call(['hisat2-build','-p','6' ,'--ss', data_dir +ss_file, '--exon', data_dir +exon_file, l_ref.split('.gz')[0], l_hisat2_index])
        """
        # Add data to database
        self.dt = self.dt.append(tobj.retDataRow(), ignore_index=True)
        writer = pandas.ExcelWriter(self.refDatabase)
        self.dt.to_excel(writer)
        writer.save()
        
    def deleteDatabase(self, databaseID):
        # Recreate database
        self.dt = pandas.read_excel(self.refDatabase)
        
        if len(self.dt[self.dt.DatabaseID == databaseID]) == 0:
            print(databaseID + ' not found. Cannot delete')
            raise KeyError
        
        self.dt = self.dt[self.dt.DatabaseID != databaseID]
        writer = pandas.ExcelWriter(self.refDatabase)
        self.dt.to_excel(writer)
        writer.save()

        # Delete directory containing data
        tobj = self.retDatabaseData(databaseID)
        shutil.rmtree(tobj.basedir)

    def retDatabaseData(self, databaseID):
        if len(self.dt[self.dt.DatabaseID == databaseID]) == 0:
            print('Not a valid database. Quitting...', file = sys.stderr)
            exit()
            
        database_data = self.dt[self.dt.DatabaseID == databaseID]
        for o_data in database_data.iterrows():
            return RefObj(databaseID, o_data[1].Species, o_data[1].BaseDir, o_data[1].GTFFlag)
        
    def retDatabaseDir(self, gv, species = 'c_elegans'):
        data_dir = self.dbDirectory + species + '/' + gv + '/'
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        return data_dir

    def retGenomeVersions(self):
        return self.dt.groupby(['DatabaseID']).groups.keys()
    
    def retGenomeVersionsBySpecies(self):
        outdata = {}
        species = self.dt.groupby(['Species']).groups.keys()
        for sp in species:
            outdata[sp] = list(self.dt[self.dt.Species == sp]['DatabaseID'])
        return outdata
    
databaseInfo = RefDatabase()

class SampleObj():

    def __init__(self, species, projectID, sampleID, libraryID, paired_flag, fqfile1, fqfile2, library_type, baseDir):
        if '.fq.gz' not in fqfile1 and '.fastq.gz' not in fqfile1:
            print(fqfile1 + ' must end in fq.gz or fastq.gz', file = sys.stderr)
            raise(NameError)
        self.species = species
        self.projectID = projectID
        self.sampleID = sampleID
        self.libraryID = str(sampleID) + '_' + str(libraryID)
        self.paired_flag = paired_flag
        self.fqfile1 = baseDir + fqfile1
        if 'processed' not in self.fqfile1:
            print(self.fqfile1)
            self.fqfile1 = self.fqfile1.replace('.fastq.', '.processed.fastq.')
            self.fqfile1 = self.fqfile1.replace('.fq.', '.processed.fq.')
            print(self.fqfile1)

        if fqfile2 == '' or fqfile2 != fqfile2:
            self.two_fq_flag = False
        else:
            if '.fq.gz' not in fqfile2 and '.fastq.gz' not in fqfile2:
                print(fqfile2 + ' must end in fq.gz or fastq.gz', file = sys.stderr)
                raise(NameError)

            self.fqfile2 = baseDir + fqfile2

            if 'processed' not in self.fqfile2:
                self.fafile2 = self.fqfile2.replace('.fastq.', '.processed.fastq.')
                self.fqfile2 = self.fqfile2.replace('.fq.', '.processed.fq.')


            self.two_fq_flag = True
        self.libraryType = library_type
        self.RG_info = '@RG\\tID:' + str(self.libraryID) + '\\tSM:' + str(self.sampleID) + '\\tPL:illumina\\tLB:' + str(self.sampleID)
        
class SampleDatabase():
    def __init__(self):
        self.sampleDirectory = installed_directory + '/Samples/'

        self.sampleDatabase = self.sampleDirectory + '/AvailableSamples.xlsx'
        self.dt = pandas.read_excel(self.sampleDatabase)

    def addSample(self, newSample):

         # 1. Reread data file to make sure it hasn't changed
        self.dt = pandas.read_excel(self.sampleDatabase)

        needed_data = ['ProjectID', 'SampleID', 'Species', 'LibraryID', 'PAIRED_FLAG', 'FastQ_1', 'FastQ_2']

        # 2. Check to make sure data is included
        for datum in needed_data:
            if datum not in newSample:
                print('Data dictionary must include ' + datum + ' key to add database....', file = sys.stderr)
                raise KeyError

        # 3. Check if SampleID is already present
        if len(self.dt[(self.dt.SampleID == newSample['SampleID']) & (self.dt.LibraryID == newSample['LibraryID'])]) != 0:
            print(newSample['SampleID'] + ' already present in the database. Cannot add sample with the same SampleID and LibraryID', file = sys.stderr)
            raise KeyError

        # 4. baseDir where data will be stored
        baseDir = self.retSampleDir(newSample['Species'], newSample['ProjectID'], newSample['SampleID'])

        # 5. Move files to appropriate directories
        if not os.path.isfile(newSample['FastQ_1']):
            raise FileNotFoundError(newSample['FastQ_1'] + ' not found')
        if newSample['FastQ_2'] != '' and newSample['FastQ_2'] == newSample['FastQ_2']:
            if not os.path.isfile(newSample['FastQ_2']):
                raise FileNotFoundError(newSample['FastQ_2'] + ' not found')
        if '.gz' not in newSample['FastQ_1']:
            call(['gzip', newSample['FastQ_1']])
            newSample['FastQ_1'] += '.gz'
        if newSample['FastQ_2'] != '' and newSample['FastQ_2'] == newSample['FastQ_2']:
            if '.gz' not in newSample['FastQ_2']:
                call(['gzip', newSample['FastQ_2']])
                newSample['FastQ_2'] += '.gz'

        newfq1 = newSample['FastQ_1'].split('/')[-1].replace('.fq.gz', '.processed.fq.gz').replace('.fastq.gz', '.processed.fastq.gz')
        newfq2 = newSample['FastQ_2'].split('/')[-1].replace('.fq.gz', '.processed.fq.gz').replace('.fastq.gz', '.processed.fastq.gz')
        
        tobj = SampleObj(newSample['Species'], newSample['ProjectID'], newSample['SampleID'], newSample['LibraryID'], newSample['PAIRED_FLAG'], newfq1, newfq2, newSample['LibraryType'], baseDir)  
        
        # 6. Process reads
        if tobj.two_fq_flag:
            call(['cutadapt', '-j', str(multiprocessing.cpu_count()), '-q', '10,10', '--pair-filter=any', '-o', tobj.fqfile1, '-p', tobj.fqfile2, newSample['FastQ_1'], newSample['FastQ_2']])
            call(['fastqc', '-t', str(multiprocessing.cpu_count()), tobj.fqfile1])
            call(['fastqc', '-t', str(multiprocessing.cpu_count()), tobj.fqfile2])
        else:
            call(['cutadapt', '-j', str(multiprocessing.cpu_count()), '-q' ,'10,10', '-o', tobj.fqfile1, newSample['FastQ_1']])
            call(['fastqc', '-t', str(multiprocessing.cpu_count()), tobj.fqfile1])

        # 7. Update database
        # Add data to database
        data = {'Species': tobj.species, 'ProjectID': tobj.projectID, 'ProjectDescription': newSample['ProjectDescription'], 'SampleID': tobj.sampleID, 'SampleDescription': newSample['SampleDescription'], 'LibraryType': tobj.libraryType, 'LibrarySelection': newSample['LibrarySelection'], 'LibraryProtocol': newSample['LibraryProtocol'], 'PAIRED_FLAG': tobj.paired_flag, 'FastQ_1': tobj.originalfqfile1.replace(baseDir,''), 'FastQ_2': tobj.originalfqfile2.replace(baseDir,''), 'InstrumentModel': newSample['InstrumentModel'], 'LibraryID': tobj.libraryID, 'RunID': tobj.libraryID}
        self.dt = self.dt.append(data, ignore_index=True)
        writer = pandas.ExcelWriter(self.sampleDatabase)
        self.dt.to_excel(writer)
        writer.save()
 
    def retSampleDir(self, species, projectID, sampleID):
        data_dir = self.sampleDirectory + species + '/' + projectID + '/' + sampleID + '/'
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        return data_dir

    def retProjectData(self, projectID):
        if len(self.dt[self.dt.projectID == projectID]) == 0:
            print('Cant find project ID: ' + projectID + ' in the database....', file = sys.stderr)
            raise KeyError
            
        database_data = self.dt[self.dt.ProjectID == projectID]
        out_data = []
        for o_data in database_data.iterrows():
            baseDir = self.retSampleDir(o_data[1].Species, o_data[1].ProjectID, o_data[1].SampleID)
            out_data.append(SampleObj(o_data[1].Species, o_data[1].ProjectID, o_data[1].SampleID, o_data[1].LibraryID, o_data[1].PAIRED_FLAG, o_data[1].FastQ_1, o_data[1].FastQ_2, o_data[1].LibraryType, baseDir))
        return out_data
    
    def retSampleIDs(self):
        return self.dt.groupby(['SampleID']).groups.keys()

    def retProjectIDs(self):
        return self.dt.groupby(['ProjectID']).groups.keys()

    def retSampleData(self, sampleID):
        out_data = []
        if len(self.dt[self.dt.SampleID == sampleID]) == 0:
            print('Not a valid sampleID. Quitting...', file = sys.stderr)
            exit()
            
        sample_data = self.dt[self.dt.SampleID == sampleID]
        for o_data in sample_data.iterrows():
            baseDir = self.retSampleDir(o_data[1].Species, o_data[1].ProjectID, o_data[1].SampleID)
            out_data.append(SampleObj(o_data[1].Species, o_data[1].ProjectID, o_data[1].SampleID, o_data[1].LibraryID, o_data[1].PAIRED_FLAG, o_data[1].FastQ_1, o_data[1].FastQ_2, o_data[1].LibraryType, baseDir))
        return out_data
    
    def delete_project(self, project_id):
        self.dt = pandas.read_excel(self.sample_database)
        self.dt = self.dt[self.dt.ProjectID != project_id]
        #Add data to remove reads from Dropbox
        self.write_sample_file()

    def delete_run(self, run_id):
        self.dt = self.dt[self.dt.RunID != run_id]
        self.write_sample_file()


        
sampleInfo = SampleDatabase()

class FileMaker():
    def __init__(self):
        pass

    def genome_error(self, species, genome_version):
        print(species + '\t' + genome_version + ' does not exits in reference database. Options are:', file = sys.stderr)
        print(', '.join(list(databaseInfo.databases.keys())), file = sys.stderr)
        #sys.exit()

    def sample_error(self, sample):
        print(sample + ' does not exist in sample database. Available samples are:', file = sys.stderr)
        print(', '.join(set(list(sampleInfo.samples.keys()))), file = sys.stderr)
        print('Available projects are:', file = sys.stderr)
        print(', '.join(set(list(sampleInfo.groupings.keys()))), file = sys.stderr)
        #sys.exit()

    def del_dir(self, analysis_type, databaseID, sample):
        base_dir = output_directory + analysis_type + '/' + databaseID + '/' + sample + '/'
        if os.path.isdir(base_dir):
            shutil.rmtree(base_dir)
                     
    def ret_asb(self, analysis, databaseID, sample, bam_type):
        return self.ret_aligned_strain_bamfile(analysis, databaseID, sample, bam_type)

    def ret_aligned_strain_bamfile(self, analysis, databaseID, sample, bam_type):

        base_bam_dir = output_directory + analysis + '/' + databaseID + '/' + sample + '/'
        
        type_name = {}       
        type_name['all'] = base_bam_dir + sample + '.merged.bam'
        type_name['discordant'] = base_bam_dir + sample + '.merged.discordant.bam'
        type_name['clipped'] = base_bam_dir + sample + '.merged.clipped.bam'
        type_name['duplication'] = base_bam_dir + sample + '.merged.duplication.bam'
        type_name['inversion'] = base_bam_dir + sample + '.merged.inversion.bam'
        type_name['unmapped'] = base_bam_dir + sample + '.merged.unmapped.bam'
        type_name['chimeric'] = base_bam_dir + sample + '.merged.chimeric.bam'
     
        if bam_type not in type_name:
            print(bam_type + ' not a legitimate option. Valid options are:')
            print(type_name.keys())
            sys.exit()
            
        if not os.path.exists(base_bam_dir):
            os.makedirs(base_bam_dir)
            
        return type_name[bam_type]

    def ret_analysis_directory(self, species, genome_version, project):
        out_dir = l_analysis_directory + species + '/' + genome_version + '/' + project + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return out_dir
    
    def ret_qualimap_dir(self, other, database, sample):
        return output_directory + other + '/' + database + '/' + sample + '/Qualimap/'
    
    def ret_temp_file(self, extension = '.sam'):  
        if not os.path.isdir(retTempDirectory()):
            os.makedirs(retTempDirectory())
        if not os.path.isdir(retTempDirectory() + 'Alignments'):
            os.makedirs(retTempDirectory() + 'Alignments')
      
        return retTempDirectory() + 'Alignments/TempFile' + str(random.randint(1, 999999999)) + extension
 
    def ret_vcf_file(self, prefix, species, gv, discovery_samples, excluded_samples, genotype_samples):
        if type(discovery_samples) is str:
            discovery_samples = [discovery_samples]
        if type(excluded_samples) is str:
            excluded_samples = [excluded_samples]
        if type(genotype_samples) is str:
            genotype_samples = [genotype_samples]

        if not os.path.isdir(l_poly_directory):
            os.makedirs(l_poly_directory)
        if not os.path.isdir(l_poly_directory + species):
            os.makedirs(l_poly_directory+species)
        if not os.path.isdir(l_poly_directory + species + '/' + gv):
            os.makedirs(l_poly_directory+species + '/' + gv)
        discovery = ';'.join(sorted(discovery_samples))
        genotyped = ';'.join(sorted(genotype_samples))
        excluded = ';'.join(sorted(excluded_samples))
            
        return l_poly_directory + species + '/' + gv + '/' + prefix + '.Sp:' + species + '.GV:' +gv + '.DS:'+discovery  + '.ES:' + excluded + '.GS:'+genotyped + '.vcf'
            
file_maker = FileMaker()
