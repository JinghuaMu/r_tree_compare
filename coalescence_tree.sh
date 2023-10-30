#!/bin/bash
# Source path and direct path
# 你查找文件的根目录，即包含所有单个基因树的大文件夹
SRC_DIR="/home/tim/mixmodel/R_tree_comparison/MAFFT-AUTO-Untrimmed/"
# 你要存放的聚合的基因树的具体路径和文件名
DST_NAME="/home/tim/mixmodel/R_tree_comparison/treefile/Single_gene_tree/tree1.treefile"

# 指定一个你的data目录下的临时文件夹，不然没有读写权限
TEMP_DIR="/home/tim/mixmodel/temp/"
mkdir -p "$TEMP_DIR"

# Find and copy file (这里你用Single运行一次， 然后将Single_ 换成 Mix_ 再运行一次，记得换一个DST_NAME)
find "$SRC_DIR" -type f -regextype posix-extended -regex ".*/Single_.*\.treefile" -exec cp {} "$TEMP_DIR" \;

# Concatenate all treefiles into one
cat "$TEMP_DIR"*.treefile > "$DST_NAME"

# Clear the temp directory
rm "$TEMP_DIR"*.treefile
