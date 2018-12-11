::batch script

FOR /L %%A IN (1,1,100) DO (

Rscript --vanilla "Desktop\Real First Chapter\Code\12_social_coa_subsamples.R" %%A

batch_coancestry.bat

Rscript --vanilla "Desktop\Real First Chapter\Code\13_mrqap.R" %%A

)