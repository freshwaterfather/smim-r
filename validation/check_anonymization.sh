#!/usr/bin/env bash
# Word-boundary grep across git-tracked files (or all non-ignored files before the first
# commit) for raw site codes/names, release years and treatment vocabulary. Any hit blocks
# the commit. Years are matched as whole numbers (not digit runs inside other numbers), and
# lines that cite the source publications are exempt.
#   validation/check_anonymization.sh
cd "$(dirname "$0")/.."
# site codes are case-sensitive (the lowercase "wb" is a file mode in renv/activate.R)
pat='\b(EA|EB|WA|WB|GEORGE|AK2|YOKO|RINGO)\b|(?i:warm|cool|treatment|replicate|dissertation)|(?<![0-9.])202[56](?![0-9])'
cite='(Volponi et al\.? \(?2025|Volponi, S\. N\.,? .*2025|2024JG008225|year: 2025|SMIMfit_2\.0|Zenodo|Liddick, M\. \(2026\))'
if git rev-parse --is-inside-work-tree >/dev/null 2>&1 && [ -n "$(git ls-files)" ]; then
  files=$(git ls-files)
else
  files=$(find . -type f -not -path './.git/*' -not -path './private/*' -not -path './inputs/*' -not -path './renv/*')
fi
hits=$(echo "$files" | grep -v -E '^(private/|inputs/|PLAN\.md$|CLAUDE\.md$|HANDOFF\.md$)|\.(html|png|tiff|pdf)$' \
       | grep -v -E '^validation/check_anonymization\.sh$' \
       | xargs -d '\n' grep -n -P "$pat" 2>/dev/null | grep -v -P "$cite")
if [ -n "$hits" ]; then
  echo "ANONYMIZATION CHECK FAILED:"; echo "$hits"; exit 1
fi
echo "anonymization check passed ($(echo "$files" | wc -l) files)"
