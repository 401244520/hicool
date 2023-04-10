"cell barcode demultiplexing"

# 导入必要的模块
import gzip

# 定义函数，用于解析fastq文件中的barcode
def parse_barcode(barcode_file):
    barcode_dict = {}
    with gzip.open(barcode_file, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                # 读取barcode序列
                barcode = line.split(":")[-1]
                # 将barcode序列作为字典的键
                barcode_dict[barcode] = []
            elif line.startswith("+"):
                # 忽略fastq文件中的“+”行
                continue
            else:
                # 读取read序列并将其加入字典中对应的barcode键的列表中
                barcode_dict[barcode].append(line)
    return barcode_dict

# 定义函数，用于将每个单细胞的reads写入对应的文件中
def write_reads(barcode_dict, output_dir):
    for barcode, reads in barcode_dict.items():
        output_file = "{}/{}.fastq.gz".format(output_dir, barcode)
        with gzip.open(output_file, "wt") as f:
            for i in range(len(reads)):
                f.write(reads[i] + "\n")

# 使用示例
barcode_file = "barcodes.fastq.gz"  # 包含barcode序列的fastq文件
output_dir = "cell_reads"  # 输出目录

# 解析barcode序列
barcode_dict = parse_barcode(barcode_file)

# 将每个单细胞的reads写入对应的文件中
write_reads(barcode_dict, output_dir)
