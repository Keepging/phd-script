#!/bin/bash
# Protein Contour Wiki — Setup Script
# Run this once in your chosen directory: bash setup.sh

echo "Creating Protein Contour Wiki structure..."

# Raw sources (immutable)
mkdir -p raw/papers raw/data raw/code raw/meetings

# Wiki (AI-maintained)
mkdir -p wiki/proteins wiki/PTM_sites wiki/PTM_writers wiki/concepts
mkdir -p wiki/analyses wiki/findings wiki/hypotheses wiki/literature
mkdir -p wiki/meetings wiki/datasets wiki/methods wiki/figures wiki/gaps

# Attachments
mkdir -p attachments

# Initial index
cat > index.md << 'EOF'
# Protein Contour Wiki — Index

> Last updated: YYYY-MM-DD
> Total pages: 0

## proteins/
(empty)

## PTM_sites/
(empty)

## PTM_writers/
(empty)

## concepts/
(empty)

## analyses/
(empty)

## findings/
(empty)

## hypotheses/
(empty)

## literature/
(empty)

## meetings/
(empty)

## datasets/
(empty)

## methods/
(empty)

## figures/
(empty)

## gaps/
(empty)
EOF

# Initial log
cat > log.md << 'EOF'
# Protein Contour Wiki — Log

Append-only. Each entry: `## [YYYY-MM-DD] operation | description`

---

## [2026-05-11] init | Wiki created
- Directory structure initialized
- CLAUDE.md schema in place
- Ready for first ingest
EOF

# Git init
git init
echo "attachments/" >> .gitignore
echo ".obsidian/" >> .gitignore
echo ".DS_Store" >> .gitignore
git add .
git commit -m "init: Protein Contour Wiki"

echo ""
echo "✅ Done. Next steps:"
echo "1. Open Obsidian → Open folder as vault → select this directory"
echo "2. Install Obsidian Git plugin → set auto-commit to 10 min"
echo "3. Create GitHub private repo → git remote add origin <url> → git push"
echo "4. Copy your source files into raw/ (papers, CSVs, code)"
echo "5. Start CC in this directory: cd $(pwd) && claude"
echo "6. Tell CC: 'Read CLAUDE.md, then ingest raw/meetings/Meeting_Log_Protein_Contour.md'"
