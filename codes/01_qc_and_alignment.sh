fastp=/path/to/fastp
STAR=/path/to/STAR-2.7.10b/bin/Linux_x86_64_static/STAR
ref_mouse_fa=/path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa
ref_mouse_gtf=/path/to/Mus_musculus.GRCm38.93.gtf
threads=10

inpath=/path/to/data_directory

ls ${inpath}/*_1.fq.gz | while read r1;
do
	r2=${r1/_1.fq/_2.fq}
	out=01_fastp/$(basename $r1 _1.fq.gz)
	echo "$out\t$r1\n\t$r2"
	
  mkdir -p $out
  fq1=$out/01_fastp_r1.fq.gz
  fq2=$out/01_fastp_r2.fq.gz
	$fastp -i $r1 -I $r2 -o $fq1 -O $fq2 \
    --trim_poly_x --thread $threads \
    --html=$out/fastp.html --json=$out/fastp.json \
    1>$out/log.fastp 2>&1
done

$STAR --runMode genomeGenerate \
  --runThreadN $threads \
  --genomeDir GRCm38_STAR_2.7.10b \
  --genomeFastafiles $ref_mouse_fa \
  --sjdbGTFfile $ref_mouse_gtf

ls 01_fastp | while read sname;
do
  R1=01_fastp/${sname}/01_fastp_r1.fq.gz
	R2=01_fastp/${sname}/01_fastp_r2.fq.gz
	
	outdir=02_alignment/${sname}
  echo "$outdir/${sname}"
	mkdir -p $outdir
	
	$STAR --runThreadN $threads \
		--readFilesIn $R1 $R2 \
		--readFilesCommand gunzip -c \
		--genomeDir GRCm38_STAR_2.7.10b \
		--quantMode GeneCounts \
		--outSAMtype None \
		--outFileNamePrefix $outdir
done
