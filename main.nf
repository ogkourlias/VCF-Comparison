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
  path vcf_ref
  path vcf_ref_tbi

  output:
  path vcf_file, emit: vcf_file
  path vcf_input_tbi, emit: vcf_file_tbi
  path vcf_ref, emit: vcf_ref
  path vcf_ref_tbi, emit: vcf_ref_tbi
  path "${vcf_file.SimpleName}_headers.txt", emit: input_headers
  path "${vcf_ref.SimpleName}_headers.txt", emit: ref_headers

  script:
  """
  tabix -H ${vcf_file} | grep "^#CHROM" | sed 's/#//g' > ${vcf_file.SimpleName}_headers.txt
  tabix -H ${vcf_ref} | grep "^#CHROM" | sed 's/#//g' > ${vcf_ref.SimpleName}_headers.txt
  """
}

process chrComp {
  // publishDir "${params.outDir}", mode: 'move'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path vcf_input
  path vcf_input_tbi
  path vcf_ref
  path vcf_ref_tbi
  path input_headers
  path ref_headers

  output:
  path "${vcf_input.simpleName}.tsv", emit: output_tsv

  script:
  """
  vcf_compare.py -f ${vcf_input} -r ${vcf_ref} -n ${params.chunkSize} -ih ${input_headers} \
  -ch ${ref_headers} -chr ${vcf_input.simpleName} -o ${vcf_input.simpleName}.tsv
  """
}

process allStats {
  // publishDir "${params.outDir}", mode: 'move'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path vcf_input

  output:
  path "${vcf_input.simpleName}.tsv.gz", emit: output_tsv

  script:
  """
  vcf_stats.py -i ${vcf_input} -o ${vcf_input.simpleName}.tsv.gz
  """
}

process combineTsv {
  //publishDir "${params.outDir}", mode: 'copy'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 0

  input:
  path tsv_files

  output:
  path "compMerged.tsv", emit: tsv_merged

  script:
  """
  cat *.tsv | head -n 1  > compMerged.tsv
  cat *.tsv | sed '/^varID/d'  >> compMerged.tsv
  """
}

process indexVcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 12

  input:
  path vcf_file

  output:
  path "${vcf_file}.tbi", emit: tbi_file
  path "${vcf_file}", emit: vcf_file

  script:
  """
  tabix -p vcf ${vcf_file}
  """
}

process get_region_comp {
  publishDir "${params.outDir}", mode: 'move'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '4 GB'
  cpus 1
  maxRetries 12

  input:
  path tsv_file

  output:
  path "comp/${tsv_file}.gz"

  script:
  """
  mkdir comp
  get_region.py -comp \
  -i ${tsv_file} \
  -g ${params.gtf_file} \
  -c ${tsv_file.simpleName} \
  -n ${params.chunkSize} \
  -o comp/${tsv_file}
  gzip comp/${tsv_file}
  """
}

process get_region_all {
  publishDir "${params.outDir}", mode: 'move'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '4 GB'
  cpus 1
  maxRetries 12

  input:
  path tsv_file

  output:
  path "all/${tsv_file}"

  script:
  """
  mkdir all
  get_region.py -all \
  -i ${tsv_file} \
  -g ${params.gtf_file} \
  -c ${tsv_file.simpleName} \
  -n ${params.chunkSize} \
  -o all/${tsv_file}
  gzip all/${tsv_file}
  """
}

workflow {
  vcf_chr_files = Channel.fromPath("${params.vcf_input}/*.vcf.gz", type: 'file')
  indexVcf(vcf_chr_files)
  getHeaders(indexVcf.output.vcf_file, indexVcf.output.tbi_file, params.vcf_ref, params.vcf_ref_tbi)  
  chrComp(getHeaders.output.vcf_file, getHeaders.output.vcf_file_tbi, getHeaders.output.vcf_ref, getHeaders.output.vcf_ref_tbi,
   getHeaders.output.input_headers, getHeaders.output.ref_headers)
  allStats(vcf_chr_files)
  get_region_comp(chrComp.output)
  get_region_all(allStats.output)
  }