
cd "C:\Users\froug\Desktop\Real First Chapter"

::there are two programs named plink so added crazy path

::version for filtering

..\plink_1.9\plink --noweb --file "original2018" --indep 50 5 2 --hwe 0.05 midp --maf 0.01 --nonfounders --allow-extra-chr --out "original2018f"

..\plink_1.9\plink --make-bed --file "original2018" --extract "original2018f.prune.in" --allow-extra-chr --out "original2018p" 

..\plink_1.9\plink --bfile "original2018p" --recode A --allow-extra-chr --out "original2018p"

..\plink_1.9\plink --bfile "original2018p" --recode tab --out "original2018p" --allow-extra-chr 

..\plink_1.9\plink --file "original2018p" --freq --allow-extra-chr

::special filtered version for sequoia 

..\plink_1.9\plink --noweb --file "original2018" --indep 50 5 2 --hwe 0.05 midp --maf 0.4 --nonfounders --allow-extra-chr --out "sequoia4"

..\plink_1.9\plink --make-bed --file "original2018" --extract "sequoia4.prune.in" --allow-extra-chr --out "sequoia4p" 

..\plink_1.9\plink --bfile "sequoia4p" --recode A --allow-extra-chr --out "sequoia4p"