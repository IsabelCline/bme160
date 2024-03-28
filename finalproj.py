#print('hello world')
#'''Required filetypes: h5ad, mtx, rds, tsv.gz, txt.gz, h5, h5seurat, fastq?, bcl?, csv?'''
#!pip install scanpy
import scanpy as sc
import streamlit as st
import anndata as ad
import os
import tempfile
import shutil
import pandas as pd
from pathlib import Path
from scipy.io import mmwrite

import pandas as pd #needed?
#could do buttons in streamlit to select infile/outfile types, that way inputs are controlled
# st.download_button("Download file", file) 
class FileConverter:
    ''' Recieve user input for file, input file type and desired file output type. 
        Translate file to AnnData object in order for it to be able to access AnnData methods.
        Use Scanpy & AnnData methods to read and write file to specified output. '''
    supported_formats = ['h5ad', 'mtx', 'csv', 'txt', 'loom', 'tsv', 'h5', 'xslx']
    def __init__(self, tempdir, inext, infilename, fname, outformat, sep): 
        ''' Initialize the file type and name for both input and output. ''' #change this
        self.dir = tempdir
        self.inext = inext
        self.infilename = infilename
        self.tfname = fname
        #st.write(self.infilename)
        self.outformat = outformat
        self.outfile = None
        self.sep = sep
    
    #methods for reading/writing that aren't given in sc
    def txt_file_write(self, adataobj, separator): #change
        '''Write an AnnData object to a txt file.'''
        return adataobj.write_csvs(self.outfile, sep = separator) #reusing csv method with defined separator
    
    def mtx_file_write(self, data): #maybe write_10x_mtx?
        ''' Write an AnnData object to an mtx file, and the corresponding genes/barcodes to tsv files.'''
        #could build in a way to let user specify header, but seems standard?
        with tempfile.TemporaryDirectory() as tempdir:
            path= Path(self.infilename) #in or out of with block?
            with open(tempdir + '/' + path.stem + '_matrix.mtx', 'wb') as f:
                st.write('About to call mmwrite...')
                mmwrite(f, data.X.T)
                st.write('Called mmwrite')
        #         f.writelines(['%%MatrixMarket matrix coordinate integer general', '%', 'placeholder'])

        # #self.outfile.writelines(['%%MatrixMarket matrix coordinate integer general', 'placeholder']) #this seems standard, is this right?
        #         nonzero = 0
        #         for row in data.AnnData.n_obs: #cells
        #             for col in data.AnnData.n_var: #genes
        #                 val = AnnData.X[row, col]
        #                 if int(val) != 0:
        #                     nonzero += 1
        #                     f.write(f'{row} {col} {val}')
        #                     #need a sort function in here somewhere?
        #         #once have counted nonzero values, then can go back and write line 2
        #         f.seek(50)#sets position to right after header to write over placeholder
        #         f.write(f'{data.AnnData.n_obs} {data.AnnData.n_var} {nonzero}') #go back and rewrite after nonzero has been calculated
        
        
        #then create the tsv files with gene names and cell names?
            with open(tempdir + '/' + path.stem + '_genes.tsv', 'wb') as genes, open(tempdir + '/' + path.stem + '_barcodes.tsv', 'wb') as cells:
                #feature = data.var["feature_types"] if data.var["feature_types"] else "Gene Expression"
                if data.var["gene_ids"]:
                    pd.DataFrame({0: data.var["gene_ids"], 1: data.var_names, 2: "Gene Expression"}).to_csv(
                        genes,
                        sep="\t",
                        index=False,
                        header=False
                        )
                elif data.var["gene_symbols"]:
                    pd.DataFrame({0: data.var_names, 1: data.var["gene_symbols"], 2: "Gene Expression"}
                        ).to_csv(genes, sep="\t", index=False, header=False)
                
                pd.DataFrame(data.obs_names).to_csv(cells, sep="\t", index=False, header=False)
                # i = 0
                # for var in data.var_names.values:
                #     st.write(var)
                #     i += 1
                #     if i > 10:
                #         break

                #for obs in data.obs_names.values:
                    #st.write(obs) #is there sorting to be done here?
            mtxoutzip = shutil.make_archive("/tmp/" + path.stem, 'zip', tempdir)
        return mtxoutzip
        
        
        #pass #not really sure how to write this one. maybe utilize scipy?

    def convert_file(self):
        ''' Call the correct methods to read and write to files of different formats.
            Try handles the cases of files that are of an included type.
            Except block returns an error message to the user so as not to crash the code when the file is not one of the included formats. '''
        #just-in-case checks; shouldn't be necessary since should be checked in streamlit

        if self.outformat not in FileConverter.supported_formats:
            st.error(f'Output format {self.outformat} is not supported.')
            return 

        input_ext = os.path.splitext(self.infilename)[1][1:].lower()
        if input_ext not in FileConverter.supported_formats:
            st.error(f'Input format {input_ext} is not supported.')
            return
        
        try:
            #reading
            if input_ext == 'csv':
                data = sc.read_csv(self.tfname)
            #cont. for supported file types
            elif input_ext == 'txt':
                data = sc.read_text(self.tfname)
            elif input_ext == 'h5ad':
                #have to read file as binary/bytes?
                data = sc.read_h5ad(self.tfname)
                
                st.write(data)
            elif input_ext == 'mtx':
                #create temp directory
                #create temp files, copy files into temp files
                #import tempfile
                #import shutil
                # with tempfile.TemporaryDirectory() as temp_dir:
                #     st.write(temp_dir)
                #     temp_infile = tempfile.TemporaryFile(dir = temp_dir)
                #     temp_gf = tempfile.TemporaryFile(dir = temp_dir)
                #     temp_cf = tempfile.TemporaryFile(dir = temp_dir)

                #     shutil.copyfileobj(self.infile, temp_infile)
                #     shutil.copyfileobj(self.gf, temp_gf)
                #     shutil.copyfileobj(self.cf, temp_cf) 
                #     st.write(os.listdir(temp_dir))

                #     #shutil.move(temp_infile.name, temp_dir)
                #     #shutil.move(temp_gf, temp_dir)
                #     #shutil.move(temp_cf, temp_dir)
                #     #temp_dir = tempfile.mkdtemp()
                #     st.write(temp_dir)
                #given zip directory, zip.extractall into a tempdirectory
                

                

                    #find(matrix) and return everything that comes before
                    #path = os.path.join(temp_dir, self.infilename) #concatenates paths instead of joining them all under one dir
                data = sc.read_10x_mtx(self.dir) #method requires a pathlike obj, so have to create a temp one
            # if input_ext == 'mtx':
            #     data = sc.read_10x_mtx(self.infile)  # sc.read_10x_mtx #10x genomics formatted mtx file
            elif input_ext == 'h5':
                data = sc.read_hdf(self.tfname)  # both for reading hd5f files, which is preferred?
            # if input_ext == 'hdf5':
            #     data = sc.read_10x_h5(self.infile)  #read 10x genomics formatted hdf5 file 
            elif input_ext == 'loom':         #loom formated hd5f
                data = sc.read_loom(self.tfname)
            elif input_ext == 'xslx':
                data = sc.read_excel(self.tfname)
            elif input_ext == 'tsv':
                data = sc.read_csv(self.infile, delimiter = '\t')

            
            self.outfile = os.path.splitext(self.infilename)[0] + f'.{self.outformat}'
            
            #writing
            if self.outformat == 'csv':
                
                #create a temp directory and give to write_csvs
                with tempfile.TemporaryDirectory() as tempdir:
                    #st.write(tempdir)

                    path= Path(self.infilename)

                    data.write_csvs(tempdir, skip_data = False) #write annotations

                    #df = data.AnnData.to_df() #returns a pandas df
                    #csvpath = tempdir + '/' + path.stem + '_gex.csv'
                    #df.to_csv(csvpath) #write gene expression matrix


                    #st.write(os.listdir(tempdir))
                    
                    #st.write(path)
                    outfilename = shutil.make_archive("/tmp/" + path.stem, 'zip', tempdir)
                    #st.write(os.listdir('/tmp'))
                    #st.write(outfilename)
                    return outfilename
            
            elif self.outformat == 'txt':
                self.txt_file_write(data, self.sep)
                return self.outfile

            elif self.outformat == 'h5ad':
                data.write_h5ad(self.outfile)
                return self.outfile

            elif self.outformat == 'mtx':
                # to mtx outfile format
                return self.mtx_file_write(data)

            elif self.outformat == 'loom':
                data.write_loom(self.outfile)

            elif self.outformat == 'tsv':
                with tempfile.TemporaryDirectory() as tempdir:
                    #st.write(tempdir)
                    
                    data.write_csvs(tempdir, skip_data = False, sep = '\t')
                    #st.write(os.listdir(tempdir))
                    for tempd, subdirs, tempfiles in os.walk(tempdir):
                        #maybe print statements?
                        for tempf in tempfiles:
                            temppath = os.path.join(tempd, tempf)
                            root, ext = os.path.splitext(temppath)
                            #st.write(f'{temppath} {tempf} {root} {ext}')
                            if ext == '.csv':
                                os.rename(temppath, root + '.tsv') #path.stem instead of os.path.basename?
                    path= Path(self.infilename)
                    #st.write(path)
                    outfilename = shutil.make_archive("/tmp/" + path.stem, 'zip', tempdir)
                    #st.write(os.listdir('/tmp'))4
                    #st.write(outfilename)
                    #have to rename all csv files to tsv
                    
                    #os.walk
                    #gives list of files
                    # for every file:
                    #     if suffix = csv:
                    #         os.rename(csv file, path.stem + tsv)
                    return outfilename
            
            #return self.outfile
            

        except Exception as e:
            st.error(f'Error converting file: {e}')
def run(): #main() analog for st
    ''' Generate the applicaton interface, with a title and welcome message for user, along with other instructions. '''
    st.title('Single-Cell Gene Expression Data File Converter')
    st.write(
        '''Welcome to the single cell gene expression data
        file converter. Your file should store a 2D matrix with cells in the rows and genes in the columns (ie. something 
        like a .tsv file that contains only gene names is not suitable for this program). 
        Currently supported formats are: csv, tsv, txt, h5ad, and loom. 
        Excel (xslx) and hdf5 (h5) files are only supported on input.
        MTX is a work in progress. 
        Get started by uploading your file below. Download button is not operational currently.''')
        
    st.write('Larger files may take a while, depending on your internet speeds.')
        
    st.write('''If you are trying to convert an mtx file into another format, you will be prompted to upload your gene 
    and cell/barcode tsv files after uploading your mtx file.
    These should have the same prefix as the mtx file.''')
    input_file = st.file_uploader("Upload a file", type=['csv', 'txt', 'mtx', 'h5ad', 'loom', 'xslx', 'h5', 'tsv'])
    if input_file is not None:
        sep = None
        #import tempfile
        #create temp directory
        input_ext = os.path.splitext(input_file.name)[1][1:].lower()
        with tempfile.TemporaryDirectory() as tempdir:
            if input_ext == 'mtx':
                gfile = st.file_uploader('Upload your tsv file containing annotated genes corresponding to the uploaded mtx file', type = ['tsv'])
                cfile = st.file_uploader('Upload your tsv file containing cell barcodes corresponding to the uploaded mtx file', type = ['tsv'])
                if gfile is not None and cfile is not None:
                    f = open(tempdir + '/matrix.mtx', 'wb')
                    f.write(input_file.getbuffer())
                    f.close() #necessary?

                    gf = open(tempdir + '/genes.tsv', 'wb')
                    gf.write(gfile.getbuffer())

                    cf = open(tempdir + '/barcodes.tsv', 'wb')
                    cf.write(cfile.getbuffer())
            else:
                f = open(tempdir + '/' + 'tempfile.' + input_ext, 'wb')
                f.write(input_file.getbuffer())
            #with tempfile.NamedTemporaryFile(dir=tempdir, suffix = input_ext) as temp_input_file:
                #temp_input_file.write(input_file.getbuffer())
                #sep = ''

                # if os.path.splitext(input_file.name)[1][1:].lower() == 'mtx':
                #     gfile = st.file_uploader('Upload your tsv file containing annotated genes corresponding to the uploaded mtx file', type = ['tsv'])
                #     cfile = st.file_uploader('Upload your tsv file containing cell barcodes corresponding to the uploaded mtx file', type = ['tsv'])
                #     #make these into tempfiles as well?


            output_format = st.selectbox("Select output format", ['csv', 'txt', 'mtx', 'loom', 'h5ad', 'tsv'])
                
            if output_format == 'txt': #ignore for now
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
                #write in ability to name the genes/cells file if output_format == 'mtx'
                # elif output_format == 'mtx': #maybe move to after convert file?
                #     gfile = st.text_input('Name the tsv file that will contain the names of the genes. Default is [filename]_genes.', value = os.path.splitext(self.infilename)[0])
                #     cfile = st.text_input('Name the tsv file that will contain the cell barcodes. Default is [filename]_barcodes.', value = os.path.splitext(self.infilename)[0])
            if st.button("Convert File"):
                converter = FileConverter(tempdir, input_ext, input_file.name, f.name, output_format, sep)
                converted_file = converter.convert_file() #will either be a single file or a zip with multiple files
                st.write(converted_file)
                if converted_file:
                    path = Path(converted_file)
                    st.success(f'File converted successfully. Output file: {path.name}')
                        #from io import BytesIO
                        #filecontent = BytesIO(open(converter.outfile, 'rb').read())
                        #with open(converter.outfile, 'rb') as file:
                            #filecontent = file.read()
                        
                    st.download_button(
                        label=f"Download {path.name}",
                        #data=filecontent,
                        data=open(converted_file, 'rb'),
                        #data = converted_file,
                        file_name=path.name)
if __name__ == '__main__':
    run()
    # '''
    # Available read methods:
    # sc.read_10x_h5 #read 10x genomics formatted hdf5 file
    # sc.read_10x_mtx #10x genomics formatted mtx file
    # sc.read_h5ad
    # sc.read_csv #csv = comma separated values. essentially same as read_text but delimiter is ','
    # #read_csv available in pandas as well
    # sc.read_excel #if we want to build in xslx files, not sure #available in pandas
    # sc.read_hdf #.h5 (hdf5) file #available in pandas as well
    # sc.read_loom
    # sc.read_mtx
    # sc.read_text #can set delimiter value (default none), may be useful for files that are similar but just have different delimiters
    # #all available in anndata as well, are essentially just copied from anndata

    # Available write methods:
    # [anndata file object].write_csvs(filename) #write to csv file
    # [anndata file object].write_loom(filename)
    # [adata fileobj].write_h5ad(filename) #.write() does the same thing?

    # Will have to write methods for the ones not listed here. Some may be very similar to ones we already have so could reuse