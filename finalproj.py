print('hello world')
'''Required filetypes: h5ad, mtx, rds, tsv.gz, txt.gz, h5, h5seurat, fastq?, bcl?, csv?'''
import scanpy as sc
import streamlit as st
import os

import pandas as pd #needed?
#could do buttons in streamlit to select infile/outfile types, that way inputs are controlled
# st.download_button("Download file", file) 
class FileConverter:
    FileConverter.supported_formats = ['h5ad', 'mtx', 'csv', 'txt']
    def __init__(self, infile, outformat): 
        '''if self.filetypeidentifier(infile) == 'mtx' and self.filetypeidentifier(outfile) == 'h5ad':
            self.outfile = self.matrix_to_h5ad(infile, outfilename)'''
            #call correct method
        self.infile = infile
        self.outformat = outformat
        self.outfile = None
    
    #methods for reading/writing that aren't given in sc
    def txt_file_write(self, fileobj, separator):
        return fileobj.write_csvs(self.outfile, sep = separator) #reusing csv method with defined separator
    
    def mtx_file_write(self):
        #write header?
        self.outfile.write(f'{data.AnnData.n_obs} {data.AnnData.n_var}')
        for row in data.AnnData.n_obs: #cells
            for col in data.AnnData.n_var: #genes
                val = AnnData.X[row, col]
                self.outfile.write(f'{row} {col} {val}')
        yield self.outfile
        #then create the tsv files with gene names and cell names?
        return genenames.tsv, cellnames.tsv #or smth along those lines
        
        
        #pass #not really sure how to write this one. maybe utilize scipy?

    def convert_file(self):
        #just-in-case checks; shouldn't be necessary since should be checked in streamlit

        if self.outformat not in FileConverter.supported_formats:
            st.error(f'Output format {self.outformat} is not supported.')
            return 

        input_ext = os.path.splitext(self.infile)[1][1:].lower()
        if input_ext not in FileConverter.supported_formats:
            st.error(f'Output format {self.outformat} is not supported.')
            return
        
        try:
            #reading
            if input_ext == 'csv':
                data = sc.read_csv(self.infile)
            #cont. for supported file types
 
            self.outfile = os.path.splitext(self.infile)[0] + f'.{self.outformat}'
            
            #writing
            if self.outformat == 'csv':
                data.write_csvs(self.outfile) #could also use pd.to_csv here
            #cont. for supported file types
            if self.outformat == 'txt':
                sep = '\t' #default separator. if used this is essentially a .tsv file
                separator_selection = st.selectbox('Select separator for data (default is tab):', ['Comma', 'Tab', '1 space', '2 spaces', '3 spaces'])
                if separator_selection == 'Comma':
                    sep = ','
                elif separator_selection == '1 space':
                    sep = ' '
                elif separator_selection == '2 spaces':
                    sep = '  '
                elif separator_selection == '3 spaces':
                    sep = '   '
                
                self.txt_file_write(data, sep)
            
            

            st.success(f'File converted successfully. Output file: {self.outfile}')

        except Exception as e:
            st.error(f'Error converting file: {e}')
    
    '''
    Available read methods:
    sc.read_10x_h5 #read 10x genomics formatted hdf5 file
    sc.read_10x_mtx #10x genomics formatted mtx file
    sc.read_h5ad
    sc.read_csv #csv = comma separated values. essentially same as read_text but delimiter is ','
    #read_csv available in pandas as well
    sc.read_excel #if we want to build in xslx files, not sure #available in pandas
    sc.read_hdf #.h5 (hdf5) file #available in pandas as well
    sc.read_loom
    sc.read_mtx
    sc.read_text #can set delimiter value (default none), may be useful for files that are similar but just have different delimiters
    #all available in anndata as well, are essentially just copied from anndata

    Available write methods:
    [anndata file object].write_csvs(filename) #write to csv file
    [anndata file object].write_loom(filename)
    [adata fileobj].write_h5ad(filename) #.write() does the same thing?

    Will have to write methods for the ones not listed here. Some may be very similar to ones we already have so could reuse
    '''
    
    # mtx infile methods:
    '''def mtx_to_h5ad(self, mtxfile, h5adfilename):
        #with open mtxfile as mtx:
        inf = sc.read_mtx(mtxfile)     #could create an AnnData obj, caching data (optional) helps speed reading process 
        outf = inf.write(h5adfilename)
        return outf
    def mtx_to_h5(self, mtxfile, h5file):
        inf = sc.read_mtx(mtxfile)
        outf = inf.write(h5file) 
        return outf
    def mtx_to_h5seurat(self, mtxfile, h5seuratfile):
        pass
    def mtx_to_loom(self, mtxfile, loomfile):
        pass
    def mtx_to_rds(self, mtxfile, rdsfile):
        pass
    def mtx_to_tsvgz(self, mtxfile, tsvgzfile):
        pass
    def mtx_to_txtgz(self, mtxfile, txtgzfile):
        pass
    
    #h5ad infile methods:
    def h5ad_to_h5(self, h5adfile, h5file):
        
    def h5ad_to_h5seurat(self, h5adfile, h5seuratfile):
        pass
    def h5ad_to_loom(self, h5adfile, loomfile):
        inf = sc.read_h5ad(h5adfile) #inf is anndata file object
        outf = inf.write_loom(loomfile)
        return outf
    def h5ad_to_mtx(self, h5adfile, mtxfile):
        pass
    def h5ad_to_rds(self, h5adfile, rdsfile):
        pass
    def h5ad_to_tsvgz(self, h5adfile, tsvgzfile):
        pass
    def h5ad_to_txtgz(self, h5adfile, txtgzfile):
        pass

    def filetypeidentifier(self, file): #necessary?
        pass'''
#change
#look online for example files
#rds files are more general
#look at anndata

#from AnnData import read_h5ad (& return desired file type?)
# anndata.read_h5ad(filename, backed=None, *, as_sparse=(), as_sparse_fmt=<class 'scipy.sparse._csr.csr_matrix'>, chunk_size=6000)[source]
