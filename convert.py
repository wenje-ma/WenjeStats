import nbformat

# ====================
# 改这里!!!
# ====================

with open(r'c:\Users\18904\Github\WenjeStats\ED4DSE\02-notes\Ch04\4.2\4.2.ipynb', 'r', encoding='utf-8') as f:
    content = f.read()
start = content.index('{')
end = content.rindex('}') + 1
nb = nbformat.reads(content[start:end], as_version=4)
output = []
for cell in nb.cells:
    if cell.cell_type == 'markdown':
        output.append(cell.source)
    elif cell.cell_type == 'code':
        if cell.source.strip():
            
            # ====================
            # # 改这里!!!
            # # ====================

            output.append('```r')
            output.append(cell.source)
            output.append('```')

md_content = '\n\n'.join(output)

# ====================
# 改这里!!!
# ====================

with open(r'c:\Users\18904\Github\WenjeStats\ED4DSE\02-notes\Ch04\4.2\4.2.md', 'w', encoding='utf-8') as f:
    f.write(md_content)

print('转换完成!')
