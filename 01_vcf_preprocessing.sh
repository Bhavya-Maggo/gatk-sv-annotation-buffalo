#!/bin/bash

# ── INPUT ────────────────────────────────────────────────────────────────────
INPUT_VCF="${1}"          # e.g. 4998.vcf
SAMPLE="${2}"             # e.g. 4998
GATK="gatk-package-4.6.1.0-local.jar"

# Auto-create output names
PASSED_VCF="passed_${SAMPLE}.vcf"
FILTERED_VCF="quality_filtered_${SAMPLE}.vcf"
HQ_VCF="high_quality_${SAMPLE}.vcf"

mkdir -p qc_stats

echo " SAMPLE: ${SAMPLE}"
echo " INPUT:  ${INPUT_VCF}"

# =============================================================================
# STEP 1 — Count Total Variants
# =============================================================================
echo ""
echo "=== STEP 1: Total Variants ==="
TOTAL=$(grep -vc "^#" "${INPUT_VCF}")
echo "Total SVs in ${INPUT_VCF}: ${TOTAL}"

# =============================================================================
# STEP 2 — SV Type Breakdown
# =============================================================================
echo ""
echo "=== STEP 2: SV Type Breakdown ==="
grep -v "^#" "${INPUT_VCF}" \
  | cut -f 8 \
  | grep -o "SVTYPE=[A-Z]*" \
  | sort | uniq -c \
  | sort -rn \
  | tee qc_stats/${SAMPLE}_svtype_raw.txt

# =============================================================================
# STEP 3 — Extract PASS Variants Only
# =============================================================================
echo ""
echo "=== STEP 3: Extracting PASS variants ==="
grep -E "^#|PASS" "${INPUT_VCF}" > "${PASSED_VCF}"
PASS_COUNT=$(grep -vc "^#" "${PASSED_VCF}")
echo "PASS variants: ${PASS_COUNT} / ${TOTAL}"
echo "Removed (non-PASS): $((TOTAL - PASS_COUNT))"

# =============================================================================
# STEP 4 — Separate SV Types (from PASS VCF)
# =============================================================================
echo ""
echo "=== STEP 4: Separating SV types ==="

for SVTYPE in DEL INS DUP INV BND; do
  OUT="${SAMPLE}_${SVTYPE}.vcf"
  grep -E "^#|SVTYPE=${SVTYPE}" "${PASSED_VCF}" > "${OUT}"
  COUNT=$(grep -vc "^#" "${OUT}")
  echo "${SVTYPE}: ${COUNT} variants → ${OUT}"
done

# =============================================================================
# STEP 5 — Apply Genotype Filters (DR/DV depth)
# =============================================================================
echo ""
echo "=== STEP 5: Applying genotype depth filters ==="

if [ ! -f "${GATK}" ]; then
  echo "WARNING: GATK jar not found: ${GATK}"
  echo "Skipping genotype filter step — copy passed VCF as quality filtered"
  cp "${PASSED_VCF}" "${FILTERED_VCF}"
else
  java -jar "${GATK}" VariantFiltration \
    -V "${PASSED_VCF}" \
    -O "${FILTERED_VCF}" \
    --genotype-filter-expression "DR < 5" --genotype-filter-name "LowRefDepth" \
    --genotype-filter-expression "DV < 5" --genotype-filter-name "LowVarDepth"
  echo "Genotype filters applied → ${FILTERED_VCF}"
fi

# Verify filters applied
echo ""
echo "=== Verifying filters applied ==="
FILTER_COUNT=$(grep -c "LowRefDepth\|LowVarDepth" "${FILTERED_VCF}" 2>/dev/null || echo 0)
if [ "${FILTER_COUNT}" -gt 0 ]; then
  echo "✓ Filters applied: ${FILTER_COUNT} variants flagged"
else
  echo "⚠ No filtered variants found — check GATK output"
fi

# =============================================================================
# STEP 6 — Remove Failing Genotype Filters
# =============================================================================
echo ""
echo "=== STEP 6: Removing low-quality variants ==="
grep -v -E "LowRefDepth|LowVarDepth|LowDP|LowGQ" "${FILTERED_VCF}" > "${HQ_VCF}"
HQ_COUNT=$(grep -vc "^#" "${HQ_VCF}")
FILTERED_COUNT=$(grep -vc "^#" "${FILTERED_VCF}")
echo "Before removal: ${FILTERED_COUNT}"
echo "After removal:  ${HQ_COUNT}"
echo "Removed:        $((FILTERED_COUNT - HQ_COUNT))"

# Confirm no filtered genotypes remain
echo ""
echo "=== Confirming no filtered genotypes remain ==="
REMAIN=$(grep -cE "LowRefDepth|LowVarDepth|LowDP|LowGQ" "${HQ_VCF}" 2>/dev/null |tr -d '[:space:]' || echo 0)
if [ "${REMAIN}" -eq 0 ]; then
  echo "✓ All low-quality genotypes removed successfully"
else
  echo "⚠ WARNING: ${REMAIN} filtered genotypes still remain — check manually"
fi

# =============================================================================
# STEP 7 — Per Chromosome SV Count
# =============================================================================
echo ""
echo "=== STEP 7: SVs per chromosome (high quality) ==="
grep -v "^#" "${HQ_VCF}" \
  | cut -f 1 \
  | sort \
  | uniq -c \
  | sort -rn \
  | tee qc_stats/${SAMPLE}_per_chromosome.txt

# =============================================================================
# STEP 8 — SV Type Breakdown After Filtering
# =============================================================================
echo ""
echo "=== STEP 8: SV Type Breakdown — High Quality ==="
grep -v "^#" "${HQ_VCF}" \
  | cut -f 8 \
  | grep -o "SVTYPE=[A-Z]*" \
  | sort | uniq -c \
  | sort -rn \
  | tee qc_stats/${SAMPLE}_svtype_hq.txt

# =============================================================================
# STEP 9 — Separate HQ SV Types
# =============================================================================
echo ""
echo "=== STEP 9: Separating high-quality SV types ==="
for SVTYPE in DEL INS DUP INV BND; do
  OUT="${SAMPLE}_HQ_${SVTYPE}.vcf"
  grep -E "^#|SVTYPE=${SVTYPE}" "${HQ_VCF}" > "${OUT}"
  COUNT=$(grep -vc "^#" "${OUT}")
  echo "${SVTYPE}: ${COUNT} → ${OUT}"
done

# =============================================================================
# STEP 10 — Final QC Summary
# =============================================================================
echo ""
echo "=============================================="
echo " FINAL SUMMARY — ${SAMPLE}"
echo "=============================================="
echo " Raw total:       ${TOTAL}"
echo " PASS only:       ${PASS_COUNT}"
echo " High quality:    ${HQ_COUNT}"
echo " Removed total:   $((TOTAL - HQ_COUNT))"
echo " Output:          ${HQ_VCF}"
echo " QC stats:        qc_stats/"
echo "=============================================="

# Save summary
cat > qc_stats/${SAMPLE}_summary.txt << EOF
Sample:         ${SAMPLE}
Raw total:      ${TOTAL}
PASS only:      ${PASS_COUNT}
High quality:   ${HQ_COUNT}
Removed:        $((TOTAL - HQ_COUNT))
Output:         ${HQ_VCF}
EOF

echo "Done."