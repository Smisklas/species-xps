import re
import numpy as np
from datetime import datetime, timedelta



class DataSet():
    def __init__(self, files):
        self.files = files
        self.groups = []
        if type(files) == list:
            for file in files:
                self.load_data(file)
        else:
            self.load_data(files)

    def load_data(self, files):
        keys_dict = {0 : 'Acquisition Date',
                         1: 'Cycle',
                         2: 'Analysis Method',
                         3: 'Analyser',
                         4: 'Analyser Slit',
                         5: 'Analyser Lens',
                         6: 'Scan Mode',
                         7: 'Curves/Scan',
                         8: 'Values/Curve',
                         9: 'Dwell Time:',
                         10: 'Excitation Energy',
                         11: 'Binding Energy',
                         12: 'Pass Energy',
                         13: 'Bias Voltage',
                         14: 'Detector Voltage',
                         15: 'Eff. Workfunction',
                         16: 'Source',
                         17: 'Comment',
                         18: 'OrdinateRange',
                         19: 'Number of Scans',
                         20: 'Scan'}

        with open(files, 'r') as f:
            file = f.read()
        data_pat = re.compile(r'(?<=#\n)[\d\.\s]+(?=\n|(?:\s*\Z))', re.DOTALL)
        data = data_pat.findall(file)
        pat = re.compile(r"""[#\s]+Group:\s+(.*)\n|
                         [#\s]+Region:\s+(.*)\n|
                         [#\s]+Acquisition\s+Date:\s+(.*)\n|
                         [#\s]+Cycle:\s+(\d*)\n|
                         [#\s]+Analysis\s+Method:\s+(.*)\n|
                         [#\s]+Analyser:\s+(.*)\n|
                         [#\s]+Analyser\s+Slit:\s+(.*)\n|
                         [#\s]+Analyser\s+Lens:\s+(.*)\n|
                         [#\s]+Scan\s+Mode:\s+(.*)\n|
                         [#\s]+Curves/Scan:\s+(.*)\n|
                         [#\s]+Values/Curve:\s+(.*)\n|
                         [#\s]+Dwell\s+Time:\s+(.*)\n|
                         [#\s]+Excitation\s+Energy:\s+(.*)\n|
                         [#\s]+Binding\s+Energy:\s+(.*)\n|
                         [#\s]+Pass\s+Energy:\s+(.*)\n|
                         [#\s]+Bias\s+Voltage:\s+(.*)\n|
                         [#\s]+Detector\s+Voltage:\s+(.*)\n|
                         [#\s]+Eff.\s+Workfunction:\s+(.*)\n|
                         [#\s]+Source:\s+(.*)\n|
                         [#\s]+Comment:\s+(.*)\n|
                         [#\s]+OrdinateRange:\s+(.*)\n|
                         [#\s]+Number of Scans:\s+(.*)\n|
                         .*Scan:\s+(\d+)\n
                         """,re.X)
        header = pat.findall(file)

        shifter = 0


        for index, line in enumerate(header):
            group, region, aq_date, *argv, scan = line
            if group:
                # Check if group exists in dataset, if so, add region to group.
                loaded_groups = [group.name for group in self.groups]
                if group in loaded_groups:
                    current_group = self.groups[loaded_groups.index(group)]
                else:
                    current_group = Group(group)
                    self.groups.append(current_group)
            elif region:
                current_region = Region(region)
                current_group.add_region(current_region)
                current_region.header['Region']=region
                ii = index+1
                group, region, *argv = header[ii]
                while len(region)==0 and ii<len(header):
                    group, region, *argv = header[ii]
                    info = [i for i, e in enumerate(argv) if len(e) != 0]
                    try:

                        label = keys_dict[info[0]]
                        if label in current_region.header.keys() and label != 'Scan':
                        # Aquisition date already in header, time to get the spectrum.
                            current_region.time.append(argv[info[0]])
                            data_block = data[shifter]
                            data_lines = data_block.strip('\n').split('\n')
                            spectrum = [[],[]]
                            for line in data_lines:
                                x, y = line.split()
                                spectrum[0].append(float(x))
                                spectrum[1].append(float(y))

                            spectrum = np.array(spectrum)
                            current_region.add_sweep(spectrum)
                            shifter += 1
                        else:
                            current_region.header[label]=argv[info[0]]

                    except IndexError:
                        break
                    ii+=1

            else:
                pass

    def __str__(self):

        return_str = ''
        for index, group in enumerate(self.groups):

            return_str+=str(index)+': '+group.name+'\n'
            for region_index, region in enumerate(group.regions):
                return_str+='->'+str(region_index)+': '+region.header['Region']+'\n'
        return return_str


class Group():
    """
    The group class which is a list of regions
    """
    def __init__(self, name):
        self.name = name
        self.regions = []
    def add_region(self, region):
        self.regions.append(region)

class Region():
    """
    The regions as defined by SPECS Prodigy. This data can be multidimensional.
    """
    def __init__(self, name):
        self.name = name
        self.header={}
        self.x = np.array([])
        self.time = []
        self.data = np.array([])


    def add_sweep(self, data):

        if len(self.x)==0:
            self.x = data[0]
            self.data = data[1]
        else:
            try:
                self.data=np.vstack((self.data,data[1]))
            except ValueError:
                return


    def generate_offset_axis(self):
        """
        If the region contains multidimensional data, this method generates an
        offset axis from the timestamps of each spectrum.
        """
        self.offset = []
        offset = 0
        format_string = '%m/%d/%y %H:%M:%S'
        for timestamp in self.time:
            date, time, utc = timestamp.split()
            if not self.offset:
                self.offset.append(0)
                offset = datetime.strptime(' '.join([date,time]), format_string)
            else:
                td = datetime.strptime(' '.join([date,time]), format_string)-offset
                self.offset.append(td.total_seconds())
        self.offset = np.array(self.offset)
        # UTC = datetime.strptime(date, "%a %b %d %H:%M:%S %Y")



    def __str__(self):

        return_str = ''
        for key in self.header.keys():
            return_str += key+':\t'+self.header[key]+'\n'

        return return_str


class TangoData(object):
    """
    Loads tango data, as obtained from the beamlines, into memory.
    The shape of the data depends on the shape of the files.
    """
    def __init__(self, filename):
        self.section=""
        self.device=''
        self.attribute=''
        self.time = []
        self.data = []
        self.load(filename)

    def load(self,filename):
        with open(filename, 'r') as f:
            buffer = f.read()
        buffer = buffer.split('\n')
        for line in buffer:
            if '#' in line and not self.section:
                # get name
                junk, line = line.split('//')
                address, section,type,device,attribute = ''.join(line).split('/')
                self.section = section
                self.device = device
                self.attribute = attribute
            elif '#' in line:
                #Skip the second line for now
                pass

            elif len(line) > 0:
                #parse the data
                time, *data = line.split('\t')
                timeformat = "%Y-%m-%d_%H:%M:%S.%f"

                try:
                    self.time.append(datetime.strptime(time, timeformat))
                except ValueError:
                    self.time.append(datetime.strptime(time, "%Y-%m-%d_%H:%M:%S"))
                if len(data) == 1:
                    self.data.append(float(data[0]))
                else:
                    self.data.append(float(data))

            else:
                pass
        self.data = np.array(self.data)



if __name__ == '__main__':
    data = DataSet(['../examples/data.xy'])
