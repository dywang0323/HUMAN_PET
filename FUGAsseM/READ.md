module load module load python/3.10.13-fasrc01  
conda activate fugassem_env

## Convert the GO term from HUMANn4 to the informative GO
```
python GeneFamily2informativeGO.py go-basic.obo \
  --mapping GO_to_UniRef90.tsv \
  --informative 20 \
  --namespace CC \
  --terms-only \
  --outfile Cat_Dog_Human_GO_rename_main_lab_RM_pre_infor_CC_0924.tsv
```
