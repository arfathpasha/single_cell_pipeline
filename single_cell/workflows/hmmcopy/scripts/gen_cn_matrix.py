'''
Created on Sep 8, 2015

@author: dgrewal
'''
import warnings
import numpy as np
import pandas as pd
import os


class GenerateCNMatrix(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(self, infile, output, sep, colname, sample_id, typ):
        self.sep = sep
        self.output = output
        self.column_name = colname
        self.input = infile
        self.sample_id = sample_id
        self.type = typ

    def get_file_format(self, output):
        _, ext = os.path.splitext(output)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        else:
            warnings.warn(
                "Couldn't detect output format. extension {}".format(ext))
            return "csv"

    def read_input_data(self):
        fileformat = self.get_file_format(self.input)

        if fileformat == "h5":
            if not self.type == "hmmcopy_corrected_reads":
                raise Exception(
                    "Can only accept hmmcopy reads data in hdf format")
            data = pd.HDFStore(
                self.input, 'r')["/hmmcopy_reads/{}".format(self.sample_id)]
        else:
            data = pd.read_csv(self.input)

        return data

    @staticmethod
    def replace_missing_vals(input_df, nan_val='N/A'):
        '''
        replace NaN values with nan_val
        '''
        input_df = input_df.fillna(nan_val)

        return input_df

    def write(self, input_df):
        '''
        write the dataframe to output file
        '''

        input_df.to_csv(self.output,
                        sep=self.sep,
                        index=False)

    def read_hmmcopy_corrected_read_file(self, sample_id):
        """

        """
        column_name = self.column_name

        data = self.read_input_data()

        if column_name in data.columns:
            df = data[['chr', 'start', 'end', 'width', column_name]]
        else:
            df = data[['chr', 'start', 'end', 'width']]

            df[column_name] = float('NaN')

        df = df.rename(columns={column_name: sample_id})

        return df

    def read_gcbias_file(self, sample_id):
        """
        parses the gcbias data
        """
        column_name = self.column_name

        data = open(self.input).readlines()
        skiprows = [i for i, v in enumerate(data) if v[0] == '#' or v == '\n']

        # If the file is empty (only header no data) then return 0s (dummy
        # data)
        try:
            data = pd.read_csv(self.input, sep='\t', skiprows=skiprows)
        except pd.io.common.EmptyDataError:
            warnings.warn('No data in the GCBias output')
            # If the file is empty (only header no data) then return 0s (dummy
            # data)
            data = np.array([np.arange(100), [0] * 100]).T
            data = pd.DataFrame(data, columns=['gc', sample_id])
            return data

        data = pd.DataFrame(data[column_name])

        data['gc'] = data.index

        df = data.rename(columns={'NORMALIZED_COVERAGE': sample_id})

        df = df[['gc', sample_id]]
        return df

    def main(self):
        '''
        main function
        '''
        sample_id = self.sample_id

        if self.type == 'hmmcopy_corrected_reads':
            data = self.read_hmmcopy_corrected_read_file(sample_id)
        else:
            data = self.read_gcbias_file(sample_id)
        self.write(data)
