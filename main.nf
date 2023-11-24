nextflow.enable.dsl=2

process getChrs {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path vcf_file
  path vcf_input_tbi

  output:
  path "*.txt"

  script:
  """
  CHRS=\$(tabix -l ${vcf_file})
  for CHR in \$CHRS
  do
    touch \$CHR.txt
  done
  """
}

process getHeaders {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path vcf_file
  path vcf_input_tbi
  path vcf_comp
  path vcf_comp_tbi

  output:
  path "${vcf_file.SimpleName}_headers.txt", emit: input_headers
  path "${vcf_comp.SimpleName}_headers.txt", emit: compare_headers

  script:
  """
  tabix -H ${vcf_file} | grep "^#CHROM" | sed 's/#//g' > ${vcf_file.SimpleName}_headers.txt
  tabix -H ${vcf_comp} | grep "^#CHROM" | sed 's/#//g' > ${vcf_comp.SimpleName}_headers.txt
  """
}

process chrComp {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path chr_txt
  path vcf_input
  path vcf_input_tbi
  path vcf_comp
  path vcf_comp_tbi
  path input_headers
  path compare_headers

  output:
  path "${chr_txt.simpleName}_1.csv", emit: input_csv
  path "${chr_txt.simpleName}_2.csv", emit: comp_csv

  script:
  """
  vcf_parse.py -f ${vcf_input} -c ${vcf_comp} -ih ${input_headers} -ch ${compare_headers} -chr ${chr_txt.simpleName}
  """
}

process combineCsv {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path csv_files

  output:
  path "input.csv"
  path "comp.csv"

  script:
  """
  cat *_1.csv | head -n 1  > input.csv
  cat *_1.csv | sed '/^CHROM/d'  >> input.csv
  cat *_2.csv | head -n 1  > comp.csv
  cat *_2.csv | sed '/^CHROM/d'  >> comp.csv
  """
}



workflow {
  getHeaders(params.vcf_input, params.vcf_input_tbi, params.vcf_comp, params.vcf_comp_tbi)
  chrs = getChrs(params.vcf_input, params.vcf_input_tbi).flatten()

  csv_files = chrComp(chrs, params.vcf_input, params.vcf_input_tbi, params.vcf_comp, 
    params.vcf_comp_tbi ,getHeaders.output.input_headers, getHeaders.output.compare_headers).collect()
  
  combineCsv(csv_files)
  
  
  }