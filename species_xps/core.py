import re
import numpy as np
from datetime import datetime, timedelta
from lmfit.models import GaussianModel

class TemperatureLog():
    def __init__(self, file):
        self.files = self.clean_filename(file)
        self.temps = []
        self.time = []
        self.load_data(file)

    def clean_filename(self, filename):
        if "/" in filename:
            #it is a partial parth
            file = filename.split('/')[-1]
        else:
            file = filename

        return file

    def load_data(self, file):
        #Loads a SPECIES temperature log and parses the data
        with open(file,'r') as f:
            data = f.read()
        data = data.strip('\n').split('\n')
        self.time = [line.split('\t')[0] for line in data]
        self.temps = np.array([float(line.split('\t')[1]) for line in data])


        #Parse the times into datetime
        format_string = '%H:%M:%S'
        date_string = self.files.split('_')[-1].replace('.txt',"")
        date = datetime.strptime(date_string, '%m-%d-%Y')
        self.time = [date + (datetime.strptime(time, format_string)-date) for time in self.time if len(time)>0]

    def match_to_data(self, to_match):
        """
        This function matches the timestamps of the temperatures and returns them.
        """
        point_to_get = []
        for point in to_match:
            shifter=[]
            for item in self.time:
                shifter.append((point - item).seconds)
            array = np.asarray(shifter)
            point_to_get.append(self.temps[np.abs(array).argmin()])
        return point_to_get



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


    def generate_time_axis(self):
        """
        Generates an axis of datetime objects from the timestamps of
        multidimensional data.
        """
        format_string = '%m/%d/%y %H:%M:%S'
        tmp_time = []
        for timestamp in self.time:
            date, time, utc = timestamp.split()
            tmp_time.append(datetime.strptime(' '.join([date,time]), format_string))

        self.time = tmp_time

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

class Align():
    """
    This is a class that aligns and normalises its constituents along a chosen energy region (default e_center+/-1 eV)
    """
    def __init__(self, datalist, e_center, minmax= (1,-1)):
        self.datalist = datalist
        self.e_center = e_center
        self.minmax = minmax
        self.stats = []
    def align(self):
        #Aligns and normalises the data to the specified region
        for data in self.datalist:
            maximum = max(data.data[(data.x > self.e_center + self.minmax[1]) & (data.x < self.e_center + self.minmax[0])])
            center = self.fit_gaussian(data, maximum)
            self.stats.append((center, maximum))

    def fit_gaussian(self, data, maximum):
        #Fits a gaussian model to a subset of the data to be aligned
        gauss = GaussianModel()
        pars = gauss.make_params(center = self.e_center, amplitude = maximum, sigma = 1)
        out = gauss.fit(data.data[(data.x > self.e_center +self.minmax[1]) & (data.x < self.e_center + self.minmax[0])]
                        , pars, x=data.x[(data.x > self.e_center +self.minmax[1]) & (data.x < self.e_center + self.minmax[0])])
        return out.params['center'].value

    def plot(self):
        if not self.stats:
            print('Not aligned, stopping...')
            return
        for stat, data in zip(self.stats, self.datalist):
            e_loc, maximum = stat
            plt.plot(data.x+(self.stats[0][0]-e_loc), data.data/maximum, label = data.name)

        plt.legend()

class ExportedIgorData(object):
    """
    test class for the exported igor data. To be included in the species-xps routine.
    """
    def __init__(self, filename, headerfile=None):
        self.file=filename.replace('exported/','').replace('.txt','')
        self.x = np.array([])
        self.data = np.array([])
        self.background = []
        self.data_wb = []
        self.weight=[]
        self.fit = None
        self.x_cal = None
        self.header = dict()
        self.load(filename)
        if headerfile:
            self.read_header(headerfile)


    def load(self, filename):
        data = np.genfromtxt(filename, delimiter=',',skip_header=1)
        data = np.transpose(data)
        self.x = data[0]
        self.data = data[1:]


    def collapse(self):
        self.data = np.sum(self.data, axis=0)


    def baseline(self, limits):
        #Removes a polynomial background from the data
        self.weight = np.array(self.x)
        self.weight[:] = 1
        if type(limits) == list:
            for region in limits:
                self.weight[(self.x < max(region)) & (self.x > min(region))] = 0
        else:
            high, low = limits
            self.weight[(self.x < max(limits)) & (self.x > min(limits))] = 0

        fit = np.polyfit(self.x, self.data, 4, w=self.weight)
        self.background = np.poly1d(fit)
        self.data_wb = self.data-self.background(self.x)


    def calibate(self, shift):
        # shifts the x scale according to shift
        self.x_cal = self.x+shift


    def fit_mod(self, mod, pars):
        # Fits a lmfit model to the data contained.
        if self.data_wb and self.x_cal:
            self.fit = mod.fit(self.data_wb, pars, x=self.x_cal)
        elif self.data_wb and not self.x_cal:
            self.fit = mod.fit(self.data_wb, pars, x=self.x)
        elif self.x_cal and not self.data_wb:
            self.fit = mod.fit(self.data, pars, x=self.x_cal)
        else:
            self.fit = mod.fit(self.data, pars, x=self.x)


    def read_header(self, headerfile):
        # Loads the header as exported from Igor Pro
        with open(headerfile,'r') as f:
            try:
                header = f.read()
            except:
                print('Unreadable file....')
                return
        header = header.replace(r'\\','/').split(r'\r')
        for index, line in enumerate(header):
            if index < 2:
                pass
            else:
                try:
                    label, info = line.split('=')
                    self.header[label]=info
                except:
                    pass


if __name__ == '__main__':
    data = DataSet(['../examples/data.xy'])
