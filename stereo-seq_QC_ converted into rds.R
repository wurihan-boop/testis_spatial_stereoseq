##########################

conda create --name st python=3.8
conda activate st
pip install --upgrade pip
memba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults
 
mamba activate  st

python
import stereo as st
import warnings
warnings.filterwarnings('ignore')

##Read spatial GEF files with the `read_gef` function of the stereopy library.
data_path ="/home/wurihan/wrh/data/st_data_tissue.gef/tissue.gef/A03700A6.tissue.gef"

st.io.read_gef_info(data_path)
##Load data to generate a StereoExpData object.
data = st.io.read_gef(
  file_path=data_path,
  bin_type='bins',
  bin_size=80,
  is_sparse=True,
)

data.tl.cal_qc()
data.tl.raw_checkpoint()

# remember to set flavor as seurat
adata = st.io.stereo_to_anndata(data,flavor='seurat',output='/home/wurihan/wrh/data/st_data_rds/A03700A6_DT_bin80_seurat_out.h5ad')

##########h5ad to rds##############
##The output .h5ad could be converted into .rds file by h5ad2rds.R.
mamba activate st
#need to download Rscript 
https://github.com/STOmics/Stereopy/blob/work/docs/source/_static/h5ad2rds.R
Rscript /home/wurihan/wrh/script/h5ad2rds.R --infile /home/wurihan/wrh/data/st_data_rds/A03700A6_DT_bin80_seurat_out.h5ad --outfile /home/wurihan/wrh/data/st_data_rds/A03700A6_DT_bin80_seurat.rds

