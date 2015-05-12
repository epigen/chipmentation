import os
import pandas as pd

root = "/media/afr/cemm-backup/chipmentation/annotation/"
annot = os.path.join(root, "hg19.cage_peak_coord_robust.TATA_Annotated.bed")
expr = os.path.join(root, "hg19.cage_peak_K562_normexpression.tsv")

annot = pd.read_csv(annot, header=None, sep="\t")
annot.columns = ["chrom", "start", "end", "score", "00Annotation", "strand", "1", "2", "color", "TATA", "CPG"]

expr = pd.read_csv(expr, sep="\t")

# Merge DFs
merged = annot.merge(expr)
merged.columns = ['chrom', 'start', 'end', 'score', 'name', 'strand', '1', '2', 'color', 'TATA', 'CPG', 'description', 'association_with_transcript', 'uniprot_id', 'exprRep1', 'exprRep2', 'exprRep3']

# Filter columns out
merged = merged[['chrom', 'start', 'end', 'score', 'name', 'strand', 'TATA', 'CPG', 'exprRep1', 'exprRep2', 'exprRep3']]

merged.to_csv(os.path.join(root, "hg19.cage_peak.robust.TATA_Annotated.expr_Annotated.tsv"), sep="\t", index=None)

# Filter expressed TSSs
mExpr = merged[(sum([merged['exprRep1'], merged['exprRep2'], merged['exprRep3']]) / 3) > 3]

mExpr.to_csv(os.path.join(root, "hg19.cage_peak.robust.TATA_Annotated.expr_Annotated.K562_Expressed.tsv"), sep="\t", index=None, header=None)
