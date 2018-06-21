from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import retTempDirectory
from ftplib import FTP, error_perm
from subprocess import call, Popen, PIPE
from collections import defaultdict
import os, sys, xlrd, pandas


class DBMaker():
    def __init__(self, downloads, local_database_file, removes):
        self.speciesSH = {'c_elegans': '', 'c_briggsae': 'cb', 'c_remanei':'cr'}
        self._qualityControl(downloads, local_database_file, removes)
        for download in self.goodDownloads:
            self._downloadFiles(download)
        if self.localDatabaseFile is not None:
            self._addLocalDatabase()
        for remove in self.goodRemoves:
            print('Deleting ' + remove, file = sys.stderr)
            DI.deleteDatabase(remove)
        
    def _qualityControl(self, downloads, local_database_file, removes):
        self.goodDownloads = set()
        self.goodRemoves = set()
        self.localDatabaseFile = None
        badDownloads = set()
        badRemoves = set()
        bad_flag = False # We use this flag so that we can analyze all data and report all the problems we can find.
        
        # No data provided; print helpful information
        if downloads is None and local_database_file is None and removes is None:
            availableGVs = DI.retGenomeVersionsBySpecies()
            print('Databases already downloaded/local that are available for analysis are:')
            for sp in sorted(availableGVs.keys()):
                print('   ' + sp + ': ' + ','.join(sorted(availableGVs[sp])))

            self._identifyWormbaseSpeciesGVs()
            print('Available species and genome versions that can be downloaded from WormBase are:')
            for sp in sorted(list(self.downloadableSpeciesGVs.keys())):
                print('   ' + sp + ': ' + ','.join(sorted(self.downloadableSpeciesGVs[sp])))
            raise NotImplementedError

        if downloads is not None:
            self._identifyWormbaseSpeciesGVs()
            for download in downloads:
                try:
                    species, gv = download.split(',')
                except ValueError:
                    badDownloads.add(download)
                    print(download + ' must contain one and only one , separating species and genome version')
                    continue
                if not self._checkWormBaseSpeciesGV(species, gv):
                    badDownloads.add(download)
                else:
                    self.goodDownloads.add((species, gv))
                    
        if local_database_file is not None:
            if not os.path.isfile(local_database_file):
                print(local_database_file + ' does not exist')
                bad_flag = True
            else:
                self.localDatabaseFile = local_database_file
                
        if removes is not None:
            availableGVs = DI.retGenomeVersions()
            for remove in removes:
                if remove not in availableGVs:
                    badRemoves.add(remove)
                else:
                    self.goodRemoves.add(remove)

        if len(badDownloads) > 0:
            bad_flag = True
            print('User asked for download that is not reasonable: ' + ','.join(sorted(list(badDownloads))))
            print('Available species and genome versions that can be downloaded from WormBase are:')
            for sp in sorted(list(self.downloadableSpeciesGVs.keys())):
                print('   ' + sp + ': ' + ','.join(sorted(self.downloadableSpeciesGVs[sp])))

        if len(badRemoves) > 0:
            bad_flag = True
            print('User asked to delete database that is not present: ' + ','.join(sorted(list(badRemoves))))
            print('Databases already downloaded/local that are available for deletion:')
            for sp in sorted(availableGVs):
                print('   ' + sp + ': ' + ','.join(sorted(availableGVs[sp])))
  
        if bad_flag:
            raise NotImplementedError

    def _downloadFiles(self, download):
        temp_db = retTempDirectory() + 'WormbaseDownloads/'
        if not os.path.exists(temp_db):
            os.makedirs(temp_db)

        #Download ref_file
        #Find ref file
        #pdb.set_trace()
        species, gv = download
        databaseID = species + gv
        if species in self.speciesSH:
            databaseID.replace(species, self.speciesSH[species])

        # Get genomic DNA file
        ftp = FTP('ftp.wormbase.org')
        ftp.login()
        ftp.cwd('pub/wormbase/species/' + species + '/sequence/genomic/')
        t_files = ftp.nlst()
        for t_file in t_files:
            if gv + '.genomic.fa.gz' in t_file:
                ref_file = t_file
                break
            
        print('Downloading ' + ref_file + ' from WormBase', file = sys.stderr)
        with open(temp_db + ref_file, 'wb') as f:
            ftp.retrbinary('RETR ' + ref_file, f.write)

        #Download annotation file
        ftp.cwd('../../gff/')
        t_files = ftp.nlst()
        for t_file in t_files:
            if gv + '.annotations.gff3.gz' in t_file:
                ann_file = t_file
                break
        # Download from ftp site
        print('Downloading ' + ann_file + ' from WormBase', file = sys.stderr)
        with open(temp_db + ann_file, 'wb') as f:
            ftp.retrbinary('RETR ' + ann_file, f.write)

        # Download transcript file
        ftp.cwd('../sequence/transcripts')
        t_files = ftp.nlst()
        tr_file = None
        for t_file in t_files:
            if gv + '.mRNA_transcripts' in t_file:
                tr_file = t_file
                break
        if tr_file == None:
            for t_file in t_files:
                if gv + '.cds_transcripts' in t_file:
                    tr_file = t_file
                    break
        # Download from ftp site
        print('Downloading ' + tr_file + ' from WormBase', file = sys.stderr)
        with open(temp_db + tr_file, 'wb') as f:
            ftp.retrbinary('RETR ' + tr_file, f.write)
                   
        ftp.close()

        # Identify gene annotations from gff file
        if '.gz' in ann_file:
            call(['gunzip', temp_db + ann_file])
            ann_file = ann_file.replace('.gz', '')
        call(['grep', 'WormBase', temp_db + ann_file], stdout = open(temp_db + 'genes.gff', 'w'))
        
        
        # Update reference database
        data = {}
        data['Species'] = species
        data['DatabaseID'] = self.speciesSH[species] + gv
        data['GenomicFile'] = temp_db + ref_file
        data['AnnotationFile'] = temp_db + 'genes.gff'
        data['TranscriptFile'] = temp_db + tr_file
        print(data)
        DI.addDatabase(data)

        """
               print('Indexing and processing reference files', file = sys.stderr)
        call(['grep', 'WormBase', data_dir + data['GtfFile']], stdout = open(data_dir + 'genes.gff', 'w'))
        with open(data_dir + 'GeneKey.txt', 'w') as f, open(data_dir + 'Genes.gff') as g:
            print('WormBaseName\tCommonName', file = f)
            for line in g:
                if "\tgene\t" in line:
                    WB = line.split('Name=')[1].split(';')[0]
                    if 'locus=' in line:
                        com = line.split('locus=')[1].split(';')[0]
                    else:
                        com = line.split('sequence_name=')[1].split(';')[0]
                    print(WB + '\t' + com, file = f)
                """

    def _addLocalDatabase(self):
        dt = pandas.read_excel(self.localDatabaseFile)
        local_file_directory = self.localDatabaseFile.replace(self.localDatabaseFile.split('/')[-1], '')
        for row in dt.iterrows():
            data = {}
            for element in row[1].iteritems():
                if element[0] == 'DatabaseID':
                    data['DatabaseID'] = element[1]
                if element[0] == 'Species':
                    data['Species'] = element[1]
                if element[0] == 'GenomeFile':
                    data['GenomicFile'] = local_file_directory + element[1]
                if element[0] == 'AnnotationFile':
                    data['AnnotationFile'] = local_file_directory + element[1]
                if element[0] == 'TranscriptFile':
                    data['TranscriptFile'] = local_file_directory + element[1]
                if element[0] == 'AnalysisDirectory':
                    data['AnalysisDirectory'] = local_file_directory + element[1]
            DI.addDatabase(data)

    def _identifyWormbaseSpeciesGVs(self):
        try:
            self.downloadableSpeciesGVs
        except AttributeError:
            downloadableSpeciesGVs = defaultdict(set)
            ftp = FTP('ftp.wormbase.org')
            ftp.login()
            ftp.cwd('pub/wormbase/species')
            species = ftp.nlst()
            for sp in species:
                try:
                    ftp.cwd('/pub/wormbase/species/' + sp + '/gff')
                except error_perm:
                    downloadableSpeciesGVs[sp] = set(['None'])
                else:
                    t_files = ftp.nlst()
                    for t_file in t_files:
                        try:
                            downloadableSpeciesGVs[sp].add('WS' + t_file.split('.WS')[1].split('.')[0])
                        except IndexError:
                            continue
            ftp.quit()
            self.downloadableSpeciesGVs = downloadableSpeciesGVs

    def _checkWormBaseSpeciesGV(self, species, gv):
        self._identifyWormbaseSpeciesGVs()
        if species not in self.downloadableSpeciesGVs.keys():
            return False
        if gv not in self.downloadableSpeciesGVs[species]:
            return False
        return True
    

