#! /usr/bin/env nextflow

params.ncpu=7
bcfs = Channel.from(file(params.chunks).readLines()).splitCsv().map{return tuple("${it[0]}.${it[1]}.${it[2]}",file(it[3]),file(it[3]+".csi"))}

ref = Channel.value(file(params.ref))
ref_fai = Channel.value(file(params.ref+".fai"))

eagle_genetic_map=Channel.value(file(params.eagle_map))
shapeit_genetic_map=Channel.value(file("${params.shapeit_map}/${params.chromosome}.b38.gmap.gz"))

chromosome=Channel.value(params.chromosome)

process filter {
    cpus 4

	input:
	tuple prefix,path(bcf),path(bcf_csi) from bcfs

	output:
	tuple path("${prefix}.flt.bcf"),path("${prefix}.flt.bcf.csi") into filtered

	shell:
	'''
	bcftools view !{bcf} -f PASS,. -i 'MAC>1 && F_MISSING<0.1' --threads !{task.cpus} -Ou \
        | bcftools annotate -x ^FORMAT/GT,FORMAT/PL -Ob -o !{prefix}.flt.bcf --threads !{task.cpus}
	bcftools index !{prefix}.flt.bcf --threads !{task.cpus}
	'''
}

process concat1 {
	publishDir 'filtered-genotypes/', mode:'copy', overwrite:true
    
	input:
	path("*") from filtered.collect()
	val(prefix) from params.chromosome

	output:
	tuple path("${prefix}.bcf"),path("${prefix}.bcf.csi") into concatenated

	shell:
	'''
	bcftools concat $(ls -v *.bcf) -n -o !{prefix}.bcf
	bcftools index !{prefix}.bcf 
    '''
}

process split {
	input:
	tuple path("input.bcf"),path("input.bcf.csi") from concatenated
    val(chrom) from params.chromosome

	output:
	stdout regions
	tuple path("input.bcf"),path("input.bcf.csi") into concatenated2
    tuple val("${chrom}.no_beagle"),path("input.bcf"),path("input.bcf.csi") into straight_to_shapeit
    
	shell:
	"""
	chunker.py input.bcf
	"""
}

beagleme=regions.splitText().map{ tuple(it.trim(),it.trim().replace(":",".").replace("-",".")) }.combine(concatenated2)

process beagle {    
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus {task.attempt+4}
	clusterOptions = "-l h_vmem=6G"
    
    input:
    tuple region,prefix,bcf,bcf_csi from beagleme
    path("ref.fa.fai") from ref_fai

    output:
    tuple path("${prefix}.beagle.bcf"),path("${prefix}.beagle.bcf.csi") into beagled    

    shell:
    '''
    bcftools annotate -r !{region} -x ^FORMAT/PL !{bcf} -Ou | bcftools view -Oz -l1 -o input.vcf.gz
    /usr/bin/java -Xss2048k -jar !{params.BEAGLE4} \
     gl=input.vcf.gz out=out nthreads=!{task.cpus} modelscale=4.0 niterations=0
    tabix out.vcf.gz
    bcftools reheader -f ref.fa.fai out.vcf.gz | bcftools view -Ob -o !{prefix}.beagle.bcf
    bcftools index --threads !{task.cpus} !{prefix}.beagle.bcf
    '''
}

process concat_again {
    publishDir 'beagle-refined/', mode:'copy', overwrite:true
    cpus 4
    executor 'local'

    input:
    path('*') from beagled.collect()
    val(chrom) from params.chromosome

    output:
    tuple val("${chrom}.beagle"),path('beagle.bcf'),path('beagle.bcf.csi') into beagled2

    shell:
    '''
    bcftools concat -aD $(ls -v *.bcf) -Ob -o beagle.bcf --threads !{task.cpus}
    bcftools index beagle.bcf --threads !{task.cpus}
    '''
}

process split_and_filter {
    cpus 4
    input:
    tuple prefix,path(input),path(input_csi) from straight_to_shapeit.concat(beagled2)
    path("ref.fa") from ref
    path("ref.fa.fai") from ref_fai

    output:
    tuple prefix,path('out.bcf'),path('out.bcf.csi') into for_shapeit,for_eagle
    shell:
    '''
    bcftools norm -f ref.fa -m -any -Ou !{input} | bcftools view -i 'MAC>1' -Ob -o out.bcf --threads !{task.cpus}
    bcftools index --threads !{task.cpus} out.bcf
    '''
}

process shapeit4 {
    cpus 14
    clusterOptions = "-l h_vmem=4G"
    publishDir 'shapeit-phased/', mode:'copy', overwrite:true

    input:
    tuple prefix,path(input),path(input_csi) from for_shapeit
    val(chrom) from chromosome
    path("gmap.gz") from shapeit_genetic_map

    output:
    tuple val("shapeit"),prefix,path("${prefix}.bcf"),path("${prefix}.bcf.csi") into shapeit
    
    shell:
    '''
    shapeit4 \
        --input !{input} --output !{prefix}.bcf --thread !{task.cpus} --sequencing \
        --map gmap.gz \
        --log !{prefix}.shapeit4.log \
        --region !{chrom}
    bcftools index !{prefix}.bcf --threads !{task.cpus}
    '''
}

process eagle2 {
    cpus 14
    clusterOptions = "-l h_vmem=4G"
    publishDir 'eagle-phased/', mode:'copy', overwrite:true

    input:
    tuple prefix,path(input),path(input_csi) from for_eagle
    val(chrom) from params.chromosome
    path("map.txt.gz") from eagle_genetic_map

    output:
    tuple val("eagle"),prefix,path("${prefix}.bcf"),path("${prefix}.bcf.csi") into eagle

    shell:
    '''
    eagle --vcf !{input} \
        --geneticMapFile map.txt.gz \
        --numThreads !{task.cpus} --outPrefix !{prefix} --vcfOutFormat b 
    bcftools index !{prefix}.bcf --threads !{task.cpus}
    '''
}

process build_bref3 {
    cpus 7
    clusterOptions = "-l h_vmem=4G"    
    publishDir "${phaser}-phased/", mode:'copy', overwrite:true

    input:
    tuple phaser,prefix,path(panel),path(panel_csi) from eagle.concat(shapeit)

    output:
    path("${prefix}.bref3")

    shell:
    '''
    bcftools view !{panel} --threads !{task.cpus} -Ov \
    | java -XX:ParallelGCThreads=4 -Xmx12G -jar !{params.BREF} > !{prefix}.bref3
    '''
}
