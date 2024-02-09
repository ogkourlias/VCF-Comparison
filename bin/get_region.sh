#!/bin/bash 
/groups/umcg-biogen/tmp01/umcg-ogkourlias/vcf_comparsion/bin/get_region.py -i /groups/umcg-biogen/tmp01/umcg-ogkourlias/gtex_comparison/tsv_raw/complete/chr22.tsv \
    -g /groups/umcg-biogen/tmp01/apps/data/UMCG/STAR_index/gencode_44_2023/gencode.v44.primary_assembly.annotation.gtf \
    -c chr22 \
    -o test1.tsv \
    -n 10000