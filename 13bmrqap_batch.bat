::batch script

FOR /L %%A IN (1,1,100) DO (

START /W Rscript --vanilla "Desktop\Real First Chapter\Code\12_social_coa_subsamples.R" %%A

START /W batch_coancestry.bat

START /W Rscript --vanilla "Desktop\Real First Chapter\Code\13_mrqap.R" %%A

)