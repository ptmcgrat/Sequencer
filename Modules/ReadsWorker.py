from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import retTempDirectory
from subprocess import call
from xml.etree import ElementTree
import os, inspect, sys, pandas, pysam

class ReadWorker():
    def __init__(self, downloads, local_read_file, removes):
        self._qualityControl(downloads, local_read_file, removes)

        if len(self.goodAccensions) > 0:
            self._downloadData()

        if self.localReadFile is not None:
            self._addLocalReads()
                
    def _qualityControl(self, downloads, local_read_file, removes):
        self.goodAccensions = set()
        self.goodRemoves = set()
        self.localReadFile = None
        badRemoves = set()
        badAccessions = set()
        bad_flag = False

        if downloads is None and local_read_file is None:
            print('Must include downloads or local read file')
            raise NotImplementedError
        
        if downloads is not None:
            print('You have chosen to download data from ENA', file = sys.stderr)
            current_samples = SI.retSampleIDs()
            current_projects = SI.retProjectIDs()
            for accension in downloads:
                if accension in current_samples or accension in current_projects:
                    response = input(accension + ' already present in Sample Database. Type y or yes to redownload and delete old files.')
                    if response not in ['y','yes']:
                        print('Skipping ' + accension, file = sys.stderr)
                        continue
                    else:
                        print('Redownloading ' + accension, file = sys.stderr)
                        if accension in current_samples:
                            SI.deleteSample(accension)
                        elif accension in current_projects:
                            SI.deleteProject(accension)
                        self.goodAccensions.add(accension)
                else:
                    self.goodAccensions.add(accension)
        
        if local_read_file is not None:
            if not os.path.isfile(local_read_file):
                print(local_read_file + ' does not exist')
                bad_flag = True
            else:
                self.localReadFile = local_read_file

        if removes is not None:
            availableSamples = SI.retSampleIDs()
            for remove in removes:
                if remove not in availableSamples:
                    badRemoves.add(remove)
                else:
                    self.goodRemoves.add(remove)
                    
        if len(badRemoves) > 0:
            bad_flag = True
            print('User asked to delete samples that are not present: ' + ','.join(sorted(list(badRemoves))))
            print('Databases already downloaded/local that are available for deletion:')
            print(','.join(sorted(availableSamples)))
                
        if bad_flag:
            raise NotImplementedError

    def _addLocalReads(self):
        dt = pandas.read_excel(self.localReadFile)
        local_file_directory = self.localReadFile.replace(self.localReadFile.split('/')[-1], '')
        for row in dt.iterrows():
            data = {}
            for element in row[1].iteritems():
                if element[0] == 'Species':
                    data['Species'] = element[1]
                elif element[0] == 'FastQ_1':
                    data['FastQ_1'] = local_file_directory + element[1]
                elif element[0] == 'FastQ_2':
                    if element[1] != '' and element[1] == element[1]:
                        data['FastQ_2'] = local_file_directory + element[1]
                    else:
                        data['FastQ_2'] = ''
                else:
                    data[element[0]] = element[1]
            SI.addSample(data)
        
    def _downloadData(self):      
        for accension in self.goodAccensions:
            print('Downloading project: ' + accension, file = sys.stderr)
            call(['python3',os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Scripts/enaBrowserTools-master/python3/enaGroupGet.py', '-f', 'fastq', '-a', '-m', '-w', '-d', retTempDirectory(), accension])
            d = retTempDirectory()
            # Get all the runs that were downloaded to that project accension
            runs = [name for name in os.listdir(d + accension) if os.path.isdir(d + accension + '/' + name)]
            for run in runs:
                print('Done downloading. Syncing data for run ' + run + ' into sample database', file = sys.stderr)
                xml_file = d + accension + '/' + run + '/' + run + '.xml'
                t_data = self.parse_xml_file(xml_file)
                SI.addSample(t_data)
        
    def parse_xml_file(self, xml_file):
        e = ElementTree.parse(xml_file).getroot()
        e_iter = e.iter()
        data = {}
        for elem in e_iter:
            if elem.tag == 'PRIMARY_ID':
                if elem.text[0:3] in ['SRR', 'ERR']:
                    data['RunID'] = elem.text.rstrip()
            if elem.tag == 'TITLE':
                data['SampleDescription'] = elem.text
                if 'elegans' in elem.text:
                    data['Species'] = 'c_elegans'
                else:
                    data['Species'] = 'unknown'
            if elem.tag == 'EXPERIMENT_REF':
                data['LibraryID'] = elem.attrib['accession'].rstrip()
            if elem.tag == 'SPOT_LENGTH':
                data['ReadLength'] = elem.text.rstrip()    
            if elem.tag == 'INSTRUMENT_MODEL':
                data['InstrumentModel'] = elem.text.rstrip()
            if elem.text == 'ENA-STUDY':
                data['ProjectID'] = next(e_iter).text.rstrip()
            if elem.text == 'ENA-SAMPLE':
                data['SampleID'] = next(e_iter).text.rstrip()
            if elem.text == 'ENA-SPOT-COUNT':
                data['NumReads'] = next(e_iter).text.rstrip()
            if elem.text == 'ENA-BASE-COUNT':
                data['ReadBases'] = next(e_iter).text.rstrip()
        e = ElementTree.parse(xml_file.replace(data['RunID'] + '.xml', data['LibraryID'] + '.xml')).getroot()
        e_iter = e.iter()
        for elem in e_iter:

            if elem.tag == 'LIBRARY_NAME':
                data['LibraryName'] = elem.text.rstrip()            
            if elem.tag == 'LIBRARY_STRATEGY':
                data['LibraryType'] = elem.text.rstrip()
            if elem.tag == 'LIBRARY_SELECTION':
                data['LibrarySelection'] = elem.text.rstrip()
            if elem.tag == 'LIBRARY_LAYOUT':
                temp = list(elem)[0].tag.rstrip()
                if temp == 'PAIRED':
                    data['PAIRED_FLAG'] = True
                else:
                    data['PAIRED_FLAG'] = False
            if elem.tag == 'LIBRARY_CONSTRUCTION_PROTOCOL':
                data['LibraryProtocol'] = elem.text.rstrip()
                if 'elegans' in elem.text.rstrip():
                    data['Species'] = 'c_elegans'

        xml_directory = os.path.dirname(os.path.realpath(xml_file))
        if not data['PAIRED_FLAG']:
            data['FastQ_1'] = xml_directory + '/' + data['RunID'] + '.fastq.gz'
            data['FastQ_2'] = ''
        else:
            data['FastQ_1'] = xml_directory + '/' + data['RunID'] + '_1.fastq.gz'
            data['FastQ_2'] = xml_directory + '/' + data['RunID'] + '_2.fastq.gz'

        return(data)
