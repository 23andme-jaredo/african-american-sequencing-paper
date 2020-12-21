#! /usr/bin/env nextflow

params.ncpu=7

bcfs = Channel.from(file(params.chunks).readLines()).splitCsv().map{return tuple(it[0],"${it[0]}.${it[1]}.${it[2]}",file(it[3]),file(it[3]+".csi"))}

ref = Channel.value(file(params.ref))
ref_fai = Channel.value(file(params.ref+".fai"))

genetic_maps=Channel.from( 1..22 ).map{tuple("chr${it}",file("${params.shapeit_map}/chr${it}.b38.gmap.gz"))}

process filter {
    cpus 4

	input:
	tuple chromosome,prefix,file(bcf),file(bcf_csi) from bcfs//.view()

	output:
	tuple chromosome,file("${prefix}.flt.bcf"),file("${prefix}.flt.bcf.csi") into filtered

	shell:
	'''
	bcftools view !{bcf} -f PASS,. -i 'MAC>1 && F_MISSING<0.1' --threads !{task.cpus} -Ou \
        | bcftools annotate -x ^FORMAT/GT,FORMAT/PL -Ob -o !{prefix}.flt.bcf --threads !{task.cpus}
	bcftools index !{prefix}.flt.bcf --threads !{task.cpus}
	'''
}

process concat1 {
	publishDir 'filtered-genotypes/', mode:'copy', overwrite:true
    maxForks 4

	input:
	tuple chromosome,file("*"),file("*") from filtered.groupTuple(by:0)//.view()

	output:
    stdout regions
	tuple chromosome,file("${chromosome}.bcf"),file("${chromosome}.bcf.csi") into concatenated

	shell:
	'''
	bcftools concat $(ls -v *.bcf) -n -o !{chromosome}.bcf
	bcftools index !{chromosome}.bcf 
    chunker.py !{chromosome}.bcf
    '''
}

beagle_input=regions.splitText().map{ tuple(it.trim().tokenize(':').get(0),it.trim(),it.trim().replace(":",".").replace("-",".")) }.combine(concatenated,by:0)

process beagle {    
    errorStrategy 'retry'
    maxRetries 3
    cpus {task.attempt+4}
	clusterOptions = "-l h_vmem=6G"
    
    input:
    tuple chrom,region,prefix,path(bcf),path(bcf_csi) from beagle_input//.view()
    path("ref.fa") from ref
    path("ref.fa.fai") from ref_fai

    output:
    tuple chrom,file("${prefix}.beagle.bcf"),file("${prefix}.beagle.bcf.csi") into beagled    

    shell:
    '''
    bcftools annotate -r !{region} -x ^FORMAT/PL !{bcf} -Ou | bcftools view -Oz -l1 -o input.vcf.gz
    /usr/bin/java -Xss2048k -jar !{params.BEAGLE4} \
     gl=input.vcf.gz out=out nthreads=!{task.cpus} modelscale=4.0 niterations=0
    tabix out.vcf.gz
    bcftools reheader -f ref.fa.fai  out.vcf.gz | bcftools view -Ob -o !{prefix}.beagle.bcf
    bcftools index --threads !{task.cpus} !{prefix}.beagle.bcf
    '''
}

process concat_again {
    publishDir 'beagle-refined/', mode:'copy', overwrite:true
    cpus 4
    maxForks 4

    input:
    tuple chrom,file('*'),file('*') from beagled.groupTuple()
    path("ref.fa") from ref
    path("ref.fa.fai") from ref_fai

    output:
    stdout phasing_regions
    tuple chrom,file("${chrom}.bcf"),file("${chrom}.bcf.csi") into beagled2

    shell:
    '''
    bcftools concat -aD $(ls -v *.bcf) -Ou \
     | bcftools norm -f ref.fa -m -any -Ou \
     | bcftools view -i 'MAC>1 && ALT!="*"' -Ob -o !{chrom}.bcf --threads !{task.cpus}
    bcftools index !{chrom}.bcf --threads !{task.cpus}
  	chunker.py !{chrom}.bcf -w 10000000 -b 1000000
    '''
}

shapeit_input=phasing_regions.splitText().map{ tuple(it.trim().tokenize(':').get(0),it.trim(),it.trim().replace(":",".").replace("-",".")) }.combine(beagled2,by:0)

process shapeit4 {
    cpus 7
    clusterOptions = "-l h_vmem=4G"

    input:
    tuple chrom,region,prefix,path(input),path(input_csi),path(map) from shapeit_input.combine(genetic_maps,by:0)//.view()

    output:
    tuple chrom,file("${prefix}.bcf"),file("${prefix}.bcf.csi") into shapeit
    
    shell:
    '''
    shapeit4 \
        --input !{input} --output !{prefix}.bcf --thread !{task.cpus} --sequencing \
        --map !{map} --region !{region} \
        --log !{prefix}.shapeit4.log
    bcftools index !{prefix}.bcf --threads !{task.cpus}
    '''
}

process concat_yet_again {
    cpus 8
    clusterOptions = "-l h_vmem=4G"    
    publishDir 'phased/', mode:'copy', overwrite:true

    input:
    tuple prefix,path('*'),path('*') from shapeit.groupTuple(by:0)

    output:
    tuple prefix,path("${prefix}.bcf"),path("${prefix}.bcf.csi") into panel_vcfs
    tuple path("${prefix}.sites.bcf"),path("${prefix}.sites.bcf.csi") into site_vcfs
    
    shell:
    '''
    bcftools concat -lc $(ls -v *.bcf) -Ou --threads !{task.cpus} \
        | bcftools norm -d none -Ob -o !{prefix}.bcf 
    bcftools index --threads !{task.cpus} !{prefix}.bcf
    bcftools view -G !{prefix}.bcf -Ob -o !{prefix}.sites.bcf
    bcftools index !{prefix}.sites.bcf
    '''
}

process build_bref3 {
    cpus 8
    clusterOptions = "-l h_vmem=4G"    
    publishDir 'bref/', mode:'copy', overwrite:true

    input:
    tuple prefix,path('in.bcf'),path('in.bcf.csi') from panel_vcfs

    output:
    file("${prefix}.bref3")

    shell:    
    '''
    bcftools view in.bcf -l1 -Oz -o panel.vcf.gz --threads !{task.cpus}
    java -XX:ParallelGCThreads=4 -Xmx16G -jar !{params.BREF} panel.vcf.gz > !{prefix}.bref3
    '''
}

process build_sites {
    publishDir 'phased/', mode:'copy', overwrite:true
    input:
    path('*') from site_vcfs.collect()//.view()
    output:
    path('sites.vcf.gz')
    path('sites.vcf.gz.csi')
    shell:
    '''
    bcftools concat $(ls -v *.bcf) -Oz -o sites.vcf.gz
    bcftools index sites.vcf.gz
    '''
}
