#!/usr/bin/env nextflow

// Special thanks to https://github.com/wurmlab/flo/blob/master/Rakefile

params.previous = false
params.current = false
params.gff3 = false

// Optional params
params.chunks = 1024

// Process params

if( !(params.previous && params.current) ) {
  exit 1, "Specify assemblies as both --previous and --current"
}

/*
if ( !(params.gff3 )) {
  exit 1, "Specify GFF3 file"
}
*/

chunks = params.chunks

Channel.fromPath(params.previous)
  .set { previous_assembly }

Channel.fromPath(params.current)
  .set { current_assembly }

Channel.fromPath([params.previous, params.current])
  .collect()
  .set { genomes }

previous = file(params.previous)
current = file(params.current)


// Convert to 2bit and get info
process convertAndGetInfo {
  cpus 1

  input:
    set file(previous), file(current) from genomes

  output:
    set file("previous.2bit"), file("current.2bit") into bits2_up
    set file("previous.sizes"), file("current.sizes") into sizes_up

  """
  faToTwoBit ${previous} previous.2bit
  twoBitInfo previous.2bit stdout | sort -k2nr > previous.sizes

  faToTwoBit ${current} current.2bit
  twoBitInfo current.2bit stdout | sort -k2nr > current.sizes
  """
}

bits2_up
  .collect()
  .set { bits2 }

sizes_up
  .collect()
  .set { sizes }

process splitCurrent {
  cpus 1
  tag { "${current_assembly.baseName}" }

  input:
    file(current_assembly)

  output:
    file("chunk_*.fa") into large_chunks

  """
  faSplit sequence ${current_assembly} ${chunks} chunk_
  """
}

large_chunks.flatMap()
  .set{ large_chunks_flat }


process splitChunks {
  cpus 1

  input:
    file(large_chunk) from large_chunks_flat

  output:
    set file("${large_chunk.baseName}.5k.fa"), file("${large_chunk.baseName}.lft") into chunks5k
  """
  faSplit -oneFile size $large_chunk 5000 ${large_chunk.baseName}.5k \
    -lift=${large_chunk.baseName}.lft
  """
}

process blatChunks {
  cpus 1

  input:
    set file(fa5k), file(lft5k) from chunks5k

  output:
    set file("${fa5k.baseName}.psl"), file(lft5k) into psl5k

  """
  blat -noHead \
    -fastMap \
    -tileSize=12 \
    -minIdentity=99 \
    ${previous} \
    ${fa5k} \
    ${fa5k.baseName}.psl

  """
}

process doLiftUp {
  cpus 1

  input:
    set file(psl5k), file(lft5k) from psl5k
    set file(previous2bit), file(current2bit) from bits2

  output:
    file("${lft5k.baseName}.sorted") into sortedChains

  """
  liftUp -type=.psl -pslQ -nohead \
    ${lft5k.baseName}.psl.lifted ${lft5k} warn ${psl5k}
  axtChain -psl -linearGap=medium ${lft5k.baseName}.psl.lifted $previous2bit \
    $current2bit ${lft5k.baseName}.chn
  chainSort ${lft5k.baseName}.chn ${lft5k.baseName}.sorted
  """
}

sortedChains.collect()
  .set{ sorted_collected }

process mergeChains {
  cpus 1
  input:
    file 'sorted' from sorted_collected

  output:
    file("combined.chn.sorted") into combined_sorted

  """
  chainMergeSort sorted* | chainSplit run stdin -lump=1
  mv run/000.chain ./combined.chn.sorted
  """
}

process chainNet {
  publishDir "liftOverChain"

  input:
    file("sorted") from combined_sorted
    set file("previousSizes"), file("currentSizes") from sizes

  output:
    file("liftOver.chain")

  """
  chainNet $sorted $previousSizes $currentSizes combined.chn.sorted.net /dev/null
  netChainSubset combined.chn.sorted.net $sorted liftOver.chain
  """
}
