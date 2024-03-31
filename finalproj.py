import scanpy as sc
import streamlit as st
import anndata as ad
import os
import tempfile
import shutil
import pandas as pd
from pathlib import Path
from scipy.io import mmwrite
from zipfile import ZipFile
import pandas as pd 
import loompy
import re

def swap_file_ext(dir, old_ext, new_ext):
    for dpath, subdirs, tempfiles in os.walk(dir):
        for fname in tempfiles:
            fpath = os.path.join(dpath, fname)
            root, ext = os.path.splitext(fpath)
            #st.write(f'{temppath} {tempf} {root} {ext}')
            if ext == old_ext:
                os.rename(fpath, root + new_ext)

class FileConverter:
    ''' Recieve user input for file, input file type and desired file output type. 
        Translate file to AnnData object in order for it to be able to access AnnData methods.
        Use Scanpy & AnnData methods to read and write file to specified output. '''
    supported_formats = ['h5ad', 'mtx', 'csv', 'txt', 'loom', 'tsv', 'h5', 'xslx'] #zip?
    def __init__(self, inpath, inext, infilename, outdir, outext, sep): 
        ''' Initialize the file type and name for both input and output. ''' #change this
        self.inpath = inpath
        self.inext = inext
        self.infilename = infilename
        self.outdir = outdir
        self.outext = outext
        self.outpath = None
        self.sep = sep
    
    def make_out_subdir(self, sfx):
        in_path = Path(self.infilename)
        subdir_path = os.path.join(self.outdir, in_path.stem + sfx)
        if os.path.exists(subdir_path):
            shutil.rmtree(subdir_path, True)
        os.makedirs(subdir_path)
        return subdir_path
    
    def mtx_file_write(self, data): # maybe name write_10x_mtx?
        ''' Write an AnnData object to an mtx file, and the corresponding genes/barcodes to tsv files.'''

        mtx_outdir_path = self.make_out_subdir('')
        mtx_fstem = os.path.join(mtx_outdir_path, Path(self.infilename).stem)
        matrix_filename = mtx_fstem + '_matrix.mtx'
        genes_filename = mtx_fstem + '_genes.tsv'
        barcodes_filename = mtx_fstem + '_barcodes.tsv'

        with open(matrix_filename, 'wb') as f:
            #st.write('About to call mmwrite...')
            mmwrite(f, data.X.T)
            #st.write('Called mmwrite')
        with open(genes_filename, 'wb') as genes, open(barcodes_filename, 'wb') as cells:
            feature = data.var["feature_types"] if data.var.get("feature_types") else "Gene Expression"
            # write genes
            if data.var.get("gene_ids"):
                pd.DataFrame({0: data.var["gene_ids"], 1: data.var_names, 2: "Gene Expression"}
                    ).to_csv(genes, sep="\t", index=False, header=False)
            elif data.var.get("gene_symbols"):
                pd.DataFrame({0: data.var_names, 1: data.var["gene_symbols"], 2: "Gene Expression"}
                    ).to_csv(genes, sep="\t", index=False, header=False)
            else:
                pd.DataFrame({0: data.var_names, 1: data.var_names, 2: "Gene Expression"}
                    ).to_csv(genes, sep="\t", index=False, header=False)
            # write cells
            pd.DataFrame(data.obs_names).to_csv(cells, sep="\t", index=False, header=False)
        # this will create a zip file with the same name as the directory
        mtxoutzip = shutil.make_archive(mtx_outdir_path, 'zip', mtx_outdir_path)
        # remove directory now that it is zipped
        shutil.rmtree(mtx_outdir_path, True)
        return mtxoutzip
        
    def txt_sep_write(self, data, sep):
        #create a temp directory and give to write_csvs
        in_format = self.inext
        if in_format[0] == '.':
            in_format = in_format[1:]
        csv_outdir_path = self.make_out_subdir('_from_' + in_format + '_to_' + self.outext )
        data.write_csvs(csv_outdir_path, skip_data = False, sep = sep)
        if self.outext != 'csv':
            swap_file_ext(csv_outdir_path, '.csv', '.' + self.outext)

        outfilename = shutil.make_archive(csv_outdir_path, 'zip', csv_outdir_path)
        shutil.rmtree(csv_outdir_path, True)

        return outfilename

    def convert_file(self):
        ''' Call the correct methods to read and write to files of different formats.
            Try handles the cases of files that are of an included type.
            Except block returns an error message to the user so as not to crash the code when the file is not one of the included formats. '''
        #just-in-case checks; shouldn't be necessary since should be checked in streamlit

        if self.outext not in FileConverter.supported_formats:
            st.error(f'Output format {self.outext} is not supported.')
            return 

        input_ext = os.path.splitext(self.infilename)[1][1:].lower() #1: strips the . before the extension
        if input_ext not in FileConverter.supported_formats:
            st.error(f'Input format {input_ext} is not supported.')
            return
        
        try:
            #reading
            if input_ext == 'csv':
                data = sc.read_csv(self.inpath)
            
            elif input_ext == 'txt':
                data = sc.read_text(self.inpath)
            
            elif input_ext == 'h5ad':
                data = sc.read_h5ad(self.inpath)
                
            elif input_ext == 'mtx':
                # check that the input path is a directory
                if not os.path.isdir(self.inpath):
                    self.inpath = os.path.dirname(self.inpath)
                # find prefix
                prefix = re.split("matrix.mtx", self.infilename)[0]
                st.write(prefix)
                if prefix != '':
                    data = sc.read_10x_mtx(self.inpath, prefix = prefix)
                else:
                    data = sc.read_10x_mtx(self.inpath)
            
            elif input_ext == 'h5':
                data = sc.read_hdf(self.inpath)  # both for reading hd5f files, which is preferred?
            
            elif input_ext == 'loom':         #loom formated hd5f
                data = sc.read_loom(self.inpath)
            
            elif input_ext == 'xslx':
                data = sc.read_excel(self.inpath)
            
            elif input_ext == 'tsv':
                data = sc.read_csv(self.inpath, delimiter = '\t')

            
            out_filename_stem = os.path.splitext(self.infilename)[0]
            uniq_outpath = os.path.join(self.outdir, out_filename_stem + f'.{self.outext}')
            i = 1
            while os.path.exists(uniq_outpath): # takes care of the case where 2 input files share the same name
                uniq_outpath = os.path.join(self.outdir, out_filename_stem + f'_{i}.{self.outext}')
                i += 1
            self.outpath = uniq_outpath

            st.write(f' -> Writing converted file: "{self.outpath}" Output format: "{self.outext}"')
            #writing
            if self.outext == 'csv':
                self.txt_sep_write(data, ',')
                # #create a temp directory and give to write_csvs
                # with tempfile.TemporaryDirectory() as tempdir:
                #     #st.write(tempdir)

                #     path= Path(self.infilename)

                #     data.write_csvs(tempdir, skip_data = False) #write annotations

                #     #df = data.AnnData.to_df() #returns a pandas df
                #     #csvpath = tempdir + '/' + path.stem + '_gex.csv'
                #     #df.to_csv(csvpath) #write gene expression matrix


                #     #st.write(os.listdir(tempdir))
                    
                #     #st.write(path)
                #     outfilename = shutil.make_archive("/tmp/" + path.stem, 'zip', tempdir)
                #     #st.write(os.listdir('/tmp'))
                #     #st.write(outfilename)
                #     return outfilename
            
            elif self.outext == 'txt':
                return self.txt_sep_write(data, self.sep)

            elif self.outext == 'h5ad':
                return data.write_h5ad(self.outpath)
                # st.write(f'self.outfile is {self.outfile}')
                # st.write(os.getcwd())
                # st.write(os.listdir(os.getcwd()))
                # return self.outpath

            elif self.outext == 'mtx':
                # to mtx outfile format
                return self.mtx_file_write(data)

            elif self.outext == 'loom':
                return data.write_loom(self.outpath)

            elif self.outext == 'tsv':
                return self.txt_sep_write(data, '\t')
            #     with tempfile.TemporaryDirectory() as tempdir:
            #         #st.write(tempdir)
                    
            #         data.write_csvs(tempdir, skip_data = False, sep = '\t')
            #         #st.write(os.listdir(tempdir))
            #         for tempd, subdirs, tempfiles in os.walk(tempdir):
            #             #maybe print statements?
            #             for tempf in tempfiles:
            #                 temppath = os.path.join(tempd, tempf)
            #                 root, ext = os.path.splitext(temppath)
            #                 #st.write(f'{temppath} {tempf} {root} {ext}')
            #                 if ext == '.csv':
            #                     os.rename(temppath, root + '.tsv') #path.stem instead of os.path.basename?
            #         path= Path(self.infilename)
            #         #st.write(path)
            #         outfilename = shutil.make_archive("/tmp/" + path.stem, 'zip', tempdir)
            #         #st.write(os.listdir('/tmp'))4
            #         #st.write(outfilename)
            #         #have to rename all csv files to tsv
                    
            #         #os.walk
            #         #gives list of files
            #         # for every file:
            #         #     if suffix = csv:
            #         #         os.rename(csv file, path.stem + tsv)
            #         return outfilename
            
            # #return self.outfile
            

        except Exception as e:
            st.error(f'Error converting file: {e}')

def write_uploaded_fileobj(in_dir, file_name, ext, mode, file_obj):
    if file_name[-1] != '.' and ext[0] != '.':
        ext = '.' + ext
    fpath = os.path.join(in_dir, file_name + ext)
    f = open(fpath, mode)
    f.write(file_obj.getbuffer())
    f.close()
    return fpath

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
        Get started by uploading your file below.''')
        
    st.write('Larger files may take a while, depending on your internet speeds.')
        
    st.write('''If you are trying to convert an mtx file into another format, you will be prompted to upload your gene 
    and cell/barcode tsv files after uploading your mtx file.
    These should have the same prefix as the mtx file.''')
    input_file = st.file_uploader("Upload a file", type=['csv', 'txt', 'mtx', 'h5ad', 'loom', 'xslx', 'h5', 'tsv', 'zip'])
    if input_file is not None:
        sep = None
        input_ext = os.path.splitext(input_file.name)[1][1:].lower()
        #st.write(input_ext)
        if input_ext == 'mtx':
            gfile = st.file_uploader('Upload your tsv file containing annotated genes corresponding to the uploaded mtx file', type = ['tsv'])
            cfile = st.file_uploader('Upload your tsv file containing cell barcodes corresponding to the uploaded mtx file', type = ['tsv'])
        
        output_ext = st.selectbox("Select output format", ['csv', 'txt', 'mtx', 'loom', 'h5ad', 'tsv'])
        
        if output_ext == 'txt':
            separator_selection = st.selectbox('Select separator for data:', ['1 space', '2 spaces', '3 spaces'])
            if separator_selection == '1 space':
                sep = ' '
            elif separator_selection == '2 spaces':
                sep = '  '
            elif separator_selection == '3 spaces':
                sep = '   '
        
        if st.button("Convert File"):
            #converted_files = []
            with tempfile.TemporaryDirectory() as tempdir:
                # make sub-directory to save uploaded file(s) in
                in_dirpath = os.path.join(tempdir, "in")
                if not os.path.exists(in_dirpath):
                    os.makedirs(in_dirpath)
                
                # make sub-directory to save converted file(s) in
                out_dirpath = os.path.join(tempdir, "out")
                if not os.path.exists(out_dirpath):
                    os.makedirs(out_dirpath)
                st.write(f'BEFORE: Input Dir: {in_dirpath}')
                st.write(os.listdir(in_dirpath))
                
                st.write(f'BEFORE: Output Dir: {out_dirpath}')
                st.write(os.listdir(out_dirpath))

                if input_ext == 'mtx':
                    #gfile = st.file_uploader('Upload your tsv file containing annotated genes corresponding to the uploaded mtx file', type = ['tsv'])
                    #cfile = st.file_uploader('Upload your tsv file containing cell barcodes corresponding to the uploaded mtx file', type = ['tsv'])
                    if gfile is not None and cfile is not None:
                        write_uploaded_fileobj(in_dirpath, 'matrix', '.mtx', 'wt', input_file)
                        write_uploaded_fileobj(in_dirpath, 'genes', '.tsv', 'wt', input_file)
                        write_uploaded_fileobj(in_dirpath, 'barcodes', '.tsv', 'wt', input_file)

                        converter = FileConverter(in_dirpath, input_ext, input_file.name, out_dirpath, output_ext, sep)
                        converter.convert_file()
                elif input_ext == 'zip':
                    #st.write('input extension is zip, line 262')
                    with ZipFile(input_file) as myzip:
                        myzip.extractall(in_dirpath)
                        #walk tempdir recursively
                    #st.write(os.listdir(tempdir))
                    for root, dirs, files in os.walk(in_dirpath): #this should walk dirs recursively?
                        #st.write(f'{root} {dirs} {files}')
                        for filename in files:
                            #st.write('\n')
                            # st.write(f'filename: {filename}')
                            # st.write(f'root: {root}')
                            filepath = os.path.join(root, filename)
                            #st.write(f'filepath: {filepath}')

                            stem, fext = os.path.splitext(filepath)
                            #st.write(f'{stem}, {fext}')

                            if fext[1:] == 'tsv':
                                if (re.search('barcodes.tsv$', filepath) is not None) or (re.search('genes.tsv$', filepath) is not None):
                                    # st.write(re.search('barcodes.tsv', filepath))
                                    # st.write(re.search('genes.tsv', filepath))
                                    st.write(f'{filepath} was skipped for conversion because it is presumed to be a barcodes/genes tsv file correlated with an mtx file.')
                            elif fext[1:] in FileConverter.supported_formats:
                                try:
                                    st.write(f'{filepath} was recognized as a {fext} file. Converting {filepath}...')
                                    converter = FileConverter(filepath, fext, filename, out_dirpath, output_ext, sep) #filename or filepath?
                                    converter.convert_file()
                                    
                                    #st.write(converted_files)
                                except Exception as e:
                                    st.error(f'Error converting file: {e}. File {filename} was not able to be converted.')

                            #elif fext == 'mtx':
                                #try:
                            else:
                                st.write(f'{filepath} was skipped for conversion because it is not a supported format.')

                    #st.write(os.listdir(tempdir))

                else:
                    infile_path = write_uploaded_fileobj(in_dirpath, 'tempfile', input_ext, 'wb', input_file)
                    converter = FileConverter(infile_path, input_ext, input_file.name, out_dirpath, output_format, sep)
                    converter.convert_file()
                
                out_file_list = os.listdir(out_dirpath)
                st.write(out_file_list)
                num_out_files = len(out_file_list)

                # if num_out_files == 0:
                #     st.error('No file was able to be converted.')
                #     return
                converted_filename = None
                converted_file_path = None

                if num_out_files > 0: # just for printing success message
                    st.success(f'Successfully converted {num_out_files} files.')
                
                if num_out_files > 1:
                    #make available to download all files as zip directory
                    input_path_obj = Path(input_file.name)
                    #st.write(f'STEM of INPUT_PATH_OBJ: {input_path_obj.stem}')
                    name = input_path_obj.stem + '_converted'
                    converted_file_path = shutil.make_archive(name, 'zip', out_dirpath)
                    #st.write(f'FILENAME returned by MAKE_ARCHIVE: {converted_file_path}')
                    converted_path_obj = Path(converted_file_path)
                    converted_filename = converted_path_obj.name
                    
                    st.download_button(
                        label=f"Download {converted_filename}",
                        data=open(converted_file_path, 'rb'),
                        file_name=converted_filename)
                    
                    st.write('Or, download each file individually below:')

                if num_out_files != 0:
                    for key, out_file in enumerate(out_file_list):
                        #st.write(converted_file)
                        converted_file_path = os.path.join(out_dirpath, out_file)
                            
                        st.download_button(
                            label=f"Download {out_file}",
                            data=open(converted_file_path, 'rb'),
                            file_name=out_file,
                            key=key)
                
                else:
                    st.error('No file was able to be converted.')
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