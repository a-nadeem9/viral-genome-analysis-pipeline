# scripts/plot_tree_toytree.py

import toytree
import toyplot.svg
import toyplot.pdf
import toyplot.png

# Snakemake‐provided paths
tree_file = snakemake.input.tree
out_svg   = snakemake.output.svg
out_pdf   = snakemake.output.pdf
out_png   = snakemake.output.png

# Load and draw the tree
tree = toytree.tree(tree_file)
draw_output = tree.draw()
canvas = draw_output[0] if isinstance(draw_output, tuple) else draw_output

# Ensure a white background so the PNG is clear
canvas.style.update({"background-color": "white"})

# Render to all formats
toyplot.svg.render(canvas, out_svg)
toyplot.pdf.render(canvas, out_pdf)
toyplot.png.render(canvas, out_png)
