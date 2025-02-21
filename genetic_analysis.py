import os
import random
import subprocess
import pandas as pd
from Bio import SeqIO

# Step 1: 랜덤 DNA 염기서열 생성 함수
def generate_random_sequence(length):
    """실제 DNA 서열처럼 보이도록 A, C, G, T 랜덤 생성"""
    return ''.join(random.choices("ACGT", k=length))

# Step 2: FASTQ 품질 점수 생성 함수
def generate_quality_scores(length):
    """FASTQ 형식의 품질 점수 ('!' ~ 'I' 범위) 생성"""
    return ''.join(random.choices("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI", k=length))

# Step 3: FASTQ 파일 생성
def generate_fastq_file(filename, num_reads=500, read_length=150):
    """
    기본적인 FASTQ 형식으로 파일 생성
    - num_reads: 리드 개수 (기본 500개)
    - read_length: 리드 길이 (기본 150bp)
    """
    with open(filename, "w") as f:
        for i in range(num_reads):
            seq_id = f"@SEQ_{i}"  # 리드 ID
            sequence = generate_random_sequence(read_length)  # 랜덤 DNA 서열 생성
            quality = generate_quality_scores(read_length)  # 품질 점수 생성
            f.write(f"{seq_id}\n{sequence}\n+\n{quality}\n")  # FASTQ 형식 저장

# Step 4: FASTA (참조 유전체) 파일 생성
def generate_fasta_file(filename, num_sequences=1, seq_length=5000):
    """
    참조 유전체(FASTA) 파일 생성
    - num_sequences: 염색체 개수
    - seq_length: 유전체 길이
    """
    with open(filename, "w") as f:
        for i in range(num_sequences):
            f.write(f">chr{i+1}\n{generate_random_sequence(seq_length)}\n")

# Step 5: FASTQ 파일 분석 함수
def analyze_fastq(fastq_file):
    """
    FASTQ 파일을 읽어 기본적인 정보 확인
    - 리드 개수 및 길이 출력
    """
    num_reads = 0
    read_lengths = []

    with open(fastq_file, "r") as f:
        for record in SeqIO.parse(f, "fastq"):
            num_reads += 1
            read_lengths.append(len(record.seq))

    print("[INFO] FASTQ 파일 분석 결과:")
    print(f"  - 파일명: {fastq_file}")
    print(f"  - 총 리드 개수: {num_reads}")
    print(f"  - 평균 리드 길이: {sum(read_lengths) / len(read_lengths) if read_lengths else 0}")

# Step 6: BWA를 사용하여 FASTQ를 BAM으로 변환
def fastq_to_bam(fastq_file, reference_fasta, output_bam):
    """
    FASTQ 데이터를 참조 유전체와 정렬하여 BAM 파일 생성
    """
    if not os.path.exists(reference_fasta):
        print(f"[오류 발생] 참조 유전체 파일 '{reference_fasta}' 없음! 작업을 중단합니다.")
        return
    
    sam_file = "output.sam"
    sorted_bam = "sorted_" + output_bam

    os.system(f"bwa index {reference_fasta}")  # 참조 유전체 인덱싱
    os.system(f"bwa mem {reference_fasta} {fastq_file} > {sam_file}")  # FASTQ → SAM 변환
    os.system(f"samtools view -bS {sam_file} > {output_bam}")  # SAM → BAM 변환
    os.system(f"samtools sort {output_bam} -o {sorted_bam}")  # BAM 정렬
    os.system(f"samtools index {sorted_bam}")  # BAM 인덱싱

    print(f"[완료됨] BAM 파일 '{sorted_bam}' 생성 완료.")

# Step 7: VCF 변이 탐지 함수
def call_variants(bam_file, reference_fasta, output_vcf):
    """
    정렬된 BAM 파일에서 변이 탐지 후 VCF 파일 생성
    """
    try:
        result = subprocess.run(
            f"bcftools mpileup -f {reference_fasta} {bam_file} | bcftools call -mv -Ov -o {output_vcf}",
            shell=True, check=True, capture_output=True, text=True
        )
        print(f"[완료됨] VCF 파일 '{output_vcf}' 생성 완료.")
    except subprocess.CalledProcessError as e:
        print(f"[오류 발생] 변이 탐지 과정에서 오류 발생: {e.stderr}")

# Step 8: VCF 파일 읽기 및 변이 정보 확인
def parse_vcf(vcf_file):
    """
    VCF 파일을 읽어 변이 데이터를 요약
    """
    variants = []
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):  # 주석 줄 제외
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 5:
                chromosome, pos, _, ref, alt = fields[:5]
                variants.append({"Chromosome": chromosome, "Position": int(pos), "Ref": ref, "Alt": alt})

    df = pd.DataFrame(variants)
    return df

# [파일명 설정]
fastq_filename = "synthetic.fastq"
fasta_filename = "synthetic_reference.fasta"
bam_output = "synthetic_aligned.bam"
vcf_output = "synthetic_variants.vcf"

# [파일 생성 실행]
generate_fastq_file(fastq_filename, num_reads=500, read_length=150)  # 500개 리드 생성
generate_fasta_file(fasta_filename, num_sequences=1, seq_length=5000)  # 참조 유전체 생성

print(f"[완료됨] FASTQ 파일 '{fastq_filename}' 생성 완료.")
print(f"[완료됨] FASTA 파일 '{fasta_filename}' 생성 완료.")

# [분석 실행]
analyze_fastq(fastq_filename)  # FASTQ 데이터 분석
fastq_to_bam(fastq_filename, fasta_filename, bam_output)  # FASTQ → BAM 변환
call_variants("sorted_" + bam_output, fasta_filename, vcf_output)  # 변이 탐지

# [VCF 파일 분석]
if os.path.exists(vcf_output):
    variant_data = parse_vcf(vcf_output)
    print("[INFO] 변이 데이터 샘플 (최대 5개)")
    print(variant_data.head(5))
else:
    print(f"[오류 발생] VCF 파일 '{vcf_output}'이 생성되지 않았습니다.")
